//! Schema and columnar conversion for the Arrow search result.

use std::sync::{Arc, Mutex, OnceLock};

use arrow_array::builder::{Float64Builder, LargeStringBuilder, StringBuilder, UInt64Builder};
use arrow_array::{ArrayRef, RecordBatch};
use arrow_schema::{DataType, Field, Schema, SchemaRef};
use risearch::alignment::fingerprint_symbols;
use risearch::{HitSink, Result, SearchHit, Sequence, TargetRegistry};

pub(crate) fn search_result_schema() -> &'static SchemaRef {
    static SCHEMA: OnceLock<SchemaRef> = OnceLock::new();
    SCHEMA.get_or_init(|| {
        Arc::new(Schema::new(vec![
            Field::new("query_idx", DataType::UInt64, false),
            Field::new("query_name", DataType::Utf8, false),
            Field::new("target_idx", DataType::UInt64, false),
            Field::new("target_name", DataType::Utf8, false),
            Field::new("q_start", DataType::UInt64, false),
            Field::new("q_end", DataType::UInt64, false),
            Field::new("t_start", DataType::UInt64, false),
            Field::new("t_end", DataType::UInt64, false),
            Field::new("strand", DataType::Utf8, false),
            Field::new("energy", DataType::Float64, false),
            Field::new("alignment", DataType::LargeUtf8, true),
        ]))
    })
}

#[derive(Default)]
struct HitColumns {
    query_idx: UInt64Builder,
    query_name: StringBuilder,
    target_idx: UInt64Builder,
    target_name: StringBuilder,
    q_start: UInt64Builder,
    q_end: UInt64Builder,
    t_start: UInt64Builder,
    t_end: UInt64Builder,
    strand: StringBuilder,
    energy: Float64Builder,
    alignment: LargeStringBuilder,
}

/// Appends each query's hits straight into the Arrow columns as that query
/// finishes, so the full hit set is never materialized at once.
pub(crate) struct ArrowSink<'a> {
    queries: &'a [(String, Sequence)],
    store: &'a TargetRegistry,
    columns: Mutex<HitColumns>,
}

impl<'a> ArrowSink<'a> {
    pub(crate) fn new(queries: &'a [(String, Sequence)], store: &'a TargetRegistry) -> Self {
        Self {
            queries,
            store,
            columns: Mutex::new(HitColumns::default()),
        }
    }

    pub(crate) fn into_batch(self, schema: &SchemaRef) -> RecordBatch {
        let mut c = self.columns.into_inner().unwrap();
        let columns: Vec<ArrayRef> = vec![
            Arc::new(c.query_idx.finish()),
            Arc::new(c.query_name.finish()),
            Arc::new(c.target_idx.finish()),
            Arc::new(c.target_name.finish()),
            Arc::new(c.q_start.finish()),
            Arc::new(c.q_end.finish()),
            Arc::new(c.t_start.finish()),
            Arc::new(c.t_end.finish()),
            Arc::new(c.strand.finish()),
            Arc::new(c.energy.finish()),
            Arc::new(c.alignment.finish()),
        ];
        RecordBatch::try_new(schema.clone(), columns).expect("column count and types match schema")
    }
}

impl HitSink for ArrowSink<'_> {
    fn consume(&self, query_idx: usize, hits: Vec<SearchHit>) -> Result<()> {
        let mut strand_buf = [0u8; 4];
        // Rendered before locking, per HitSink's contract.
        let rendered: Vec<Option<String>> = hits
            .iter()
            .map(|h| {
                h.alignment
                    .as_ref()
                    .map(|a| fingerprint_symbols(a).collect())
            })
            .collect();
        let query_name = &self.queries[query_idx].0;

        let mut c = self.columns.lock().unwrap();
        for (h, alignment) in hits.iter().zip(&rendered) {
            let target_idx = h.target_index();
            c.query_idx.append_value(u64::from(h.query_idx));
            c.query_name.append_value(query_name);
            c.target_idx.append_value(u64::from(h.target_idx));
            c.target_name.append_value(self.store.get_name(target_idx));
            // The index format is 64-bit only, so usize never exceeds u64.
            c.q_start.append_value(h.q_start as u64);
            c.q_end.append_value(h.q_end as u64);
            c.t_start.append_value(h.t_start as u64);
            c.t_end.append_value(h.t_end as u64);
            c.strand
                .append_value(char::from(h.strand).encode_utf8(&mut strand_buf));
            c.energy.append_value(f64::from(h.energy));
            match alignment {
                Some(a) => c.alignment.append_value(a),
                None => c.alignment.append_null(),
            }
        }
        Ok(())
    }
}
