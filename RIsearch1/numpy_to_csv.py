import pandas as pd
import numpy as np
import os

def numpy_to_csv(matrix = np.random.rand(36, 36) , filename = "./modified_dsm.csv"):
    # Get the directory from the provided path
    dir_path = os.path.dirname(filename)

    # Construct the full path to 'risearch.c' in the same directory
    risearch_file_path = os.path.join(dir_path, 'risearch.c')

    # Assert if the file exists in the directory
    assert os.path.exists(risearch_file_path), f"Error: 'risearch.c' not found in the directory: {dir_path}"
    assert matrix.shape == (36, 36), f"Expected shape (36, 36), but got {matrix.shape}"
    df = pd.DataFrame(matrix)
    df.to_csv(filename, index=False, header=False)

if __name__ == "__main__":

    import ctypes
    import numpy as np


    # Load the shared library
    lib = ctypes.CDLL("./dsm.so")  # Use "dsm.dll" on Windows

    # Define a function to load and reshape the 4D array
    def load_and_reshape_array(var_name):
        """Load a 6x6x6x6 C array and reshape it into 36x36 NumPy array."""
        # Define the array type (4D short array)
        ArrayType = ctypes.c_short * 6 * 6 * 6 * 6

        # Get a reference to the array from the shared library
        c_array = ArrayType.in_dll(lib, var_name)

        # Convert to a NumPy array
        np_array = np.ctypeslib.as_array(c_array).reshape((6, 6, 6, 6))

        # Reshape it to 36x36
        np_reshaped = np_array.reshape((36, 36))

        return np_reshaped

    # Call the function for both variables
    dsm_su95_rev_woGU_pos = load_and_reshape_array("dsm_su95_rev_woGU_pos")
    dsm_extend = load_and_reshape_array("dsm_extend")

    numpy_to_csv(matrix = dsm_su95_rev_woGU_pos)
    numpy_to_csv(matrix = dsm_extend ,filename ="./modified_dsm_extend.csv")
