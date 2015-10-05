import java.util.Arrays;

public class HWK2_luderswd {
	public static void main(String[] strArgs) {
		/*
		 * Name: William Luders MacID: luderswd Student Number: 1403581
		 * Description: Implements linear algebra operations for matrices, such
		 * as adjoint, dot product, matrix multiplication. Returns inverse of
		 * product of some N matrices
		 */

		double[] doubleArgs = parseToDouble(strArgs); // parses arguments to doubles

		int numMatrices = (int) Math.round(doubleArgs[0]); // gets number of matrices as integer
		int finalNumRows = (int) Math.round(doubleArgs[1]); // gets row dimension of final square matrix
		int finalNumColumns = (int) Math.round(doubleArgs[numMatrices * 2]); // gets column dimensions of final square matrix

		doubleArgs = Arrays.copyOfRange(doubleArgs, 1, doubleArgs.length); // recreates array without first entry

		double[][][] matrices = new double[numMatrices][][]; // 1D array of 2D matrices, without data

		int[][] matrixSizes = new int[numMatrices][2]; // init matrixSizes, two entries for each of the matrices (dimensions)

		for (int i = 0; i < numMatrices; i++) {// iterate over matrix sizes
			matrixSizes[i][0] = (int) Math.round(doubleArgs[i * 2]); // finds row dimension of ith matrix
			matrixSizes[i][1] = (int) Math.round(doubleArgs[i * 2 + 1]); // finds column dimension of ith matrix
		}

		Boolean sizesAgree = true; // flag for whether matrix sizes allow for multiplication
		for (int i = 0; i < matrixSizes.length - 1; i++) { // iterate over matrices
			if (matrixSizes[i][1] != matrixSizes[i + 1][0]) // check that dimensions agree for multiplication
				sizesAgree = false; // sets flag to false, will be acted upon later
		}
		if (sizesAgree == false)
			System.out.println("Multiplication error"); // if dimensions don't agree, print "Multiplication error"
		else { // rest of calculations
			for (int i = 0; i < numMatrices; i++) { // iterate over matrices
				double[][] matrix = new double[matrixSizes[i][0]][matrixSizes[i][1]]; // create matrix with desired size
				matrices[i] = matrix; // populate array of with empty matrices of proper size
			}

			doubleArgs = Arrays.copyOfRange(doubleArgs, numMatrices * 2, doubleArgs.length); // recreates array without matrix size data

			for (int i = 0; i < numMatrices; i++) { // iterate over matrices
				for (int j = 0; j < matrixSizes[i][0]; j++) { // iterate over rows of matrices
					for (int k = 0; k < matrixSizes[i][1]; k++) // iterate over columns of matrices
						matrices[i][j][k] = doubleArgs[k]; // fills matrices with input values
					doubleArgs = Arrays.copyOfRange(doubleArgs, matrixSizes[i][1], doubleArgs.length); //delete matrix values that were just put into matrices
				}
			}

			double[][] matrixProduct = new double[matrixSizes[0][0]][matrixSizes[0][1]]; // initialize size for final matrix
			matrixProduct = (multiplyMatrices(matrices)); // product of all matrices

			double determinant = getDeterminant(matrixProduct); // calls function to get determinant of matrixProduct

			if (determinant == 0 || finalNumColumns != finalNumRows) // checks criteria for invertibility
				System.out.println("Matrix not invertible"); // duh
			else { // if matrix invertible
				double[][] inverseMatrix = new double[matrixProduct.length][matrixProduct.length]; // init inverse matrix with same size as original matrix
				inverseMatrix = scalarMult(getCofactor(matrixProduct), 1.0 / determinant); // divides each entry of inverse matrix by the determinant
				printMatrix(inverseMatrix); // prints result
			}
		}
	}

	public static double[] parseToDouble(String[] args) { // converts input string to doubles
		double[] doubleArgs = new double[args.length]; // create new array, same length as input array
		for (int i = 0; i < args.length; i++) // iterate over elements of new array
			doubleArgs[i] = Integer.parseInt(args[i]) * 1.0; // fill each element with input cast as double
		return doubleArgs; // returns new array as output
	}

	public static double dotProduct(double[] row, double[] column) { // performs dot product on two given arrays, in this case a ro and column
		double product = 0; // initialize dot product vbl and set to 0
		for (int i = 0; i < row.length; i++) // iterate over length of array
			product += row[i] * column[i]; // multiply values at index i of each vector and add to dot product
		return product; // return computed dot product
	}

	public static double[][] scalarMult(double[][] matrix, double scalar) { // multiplies given 2D matrix my a scalar quantity
		double[][] scaledMatrix = new double[matrix.length][matrix.length]; // declares new matrix to fill with scaled values
		for (int i = 0; i < matrix.length; i++) { // iterate over rows of matrix
			for (int j = 0; j < matrix.length; j++) // iterate over columns of matrix
				scaledMatrix[i][j] = scalar * matrix[i][j]; // multiply element of matrix by scalar
		}
		return scaledMatrix; // speaks for itself, doesn't it? Like all good code should
	}

	public static double[][] multiplyMatrices(double[][][] matrices) { // multiplies sequential matrices given in a 1D array of 2D matrices (total 3D)
		if (matrices.length == 1) { // if only 1 matrix .. 
			return matrices[0]; // return the single 2D matrix in the given 3D array
		} else { // if more than one matrix
			double[][] matrixProduct = new double[matrices[0].length][matrices[1][0].length]; // init new matrix to contain product of two matrices, with row dim of first matrix, column dim of second matrix
			for (int j = 0; j < matrices[0].length; j++) { // iterate over rows of first matrix
				for (int k = 0; k < matrices[1][0].length; k++) { // iterate over columns of second matrix
					matrixProduct[j][k] = dotProduct(getRow(matrices[0], j), getColumn(matrices[1], k)); // dot ith row of 1st matrix with jth row of 2nd matrix, add to matrix product
				}
			}
			matrices[1] = matrixProduct; // replace second entry in matrix with resultant product matrix
			return (multiplyMatrices(Arrays.copyOfRange(matrices, 1, matrices.length)));
			/*
			 * perform multiplication with resulting matrix and next matrix.
			 * first 2D matrix gets deleted as new 3D matrix is copied without
			 * index 0
			 */
		}
	}

	public static double getDeterminant(double[][] matrix) { // returns determinant of passed-in matrix
		double det = 0.0; // initialize and set determinant to 0
		if (matrix.length == 1) // if matrix is 1x1
			return matrix[0][0]; // return single entry of matrix. Base case
		if (matrix.length == 2) // if matrix is 2x2
			return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]; // return det = ad-bc. Another base case
		for (int i = 0; i < matrix.length; i++) // iterate over first row, perform cofactor expansion
			det += Math.pow(-1, i) * matrix[0][i] * getDeterminant(getMinor(matrix, 0, i)); // recursive call to find signed determinant of minor
		return det; // returns computed value of determinant
	}

	public static double[][] getCofactor(double[][] matrix) { // returns transposed cofactor matrix
		double[][] cofactor = new double[matrix.length][matrix[0].length]; // init cofactor matrix, same size as passed-in matrix
		for (int i = 0; i < cofactor.length; i++) { // iterate over rows of passed matrix
			for (int j = 0; j < cofactor.length; j++) { // iterate over columns of passed-in matrix 
				cofactor[j][i] = (getDeterminant(getMinor(matrix, i, j)) * (int) Math.pow(-1, i + j)); // 
				/*
				 * cofactor entry at (i, j) is determinant of signed minor
				 * matrix with row i and column j removed. Also transposes to
				 * find adjoint matrix, as we actually assign result to (jth,
				 * ith) entry of cofactor matrix
				 */
			}
		}
		return cofactor; // return the computed transposed cofactor matrix
	}

	public static double[][] getMinor(double[][] matrix, int rowRemoved, int columnRemoved) { // returns matrix with specified row and column removed
		double[][] minor = new double[matrix.length - 1][matrix.length - 1]; // initialize matrix with one less row and column
		int rowShift = 0; // determines if removed row has been passed
		int colShift = 0; // determines if removed column has been passed
		for (int row = 0; row < matrix.length; row++) { // iterate over rows of matrix
			colShift = 0; // resets colShift after each pass over row
			if (row == rowRemoved) { // when loop comes to desired row to be removed
				rowShift = 1; // don't do anything with this row, and indicate row has been reached
			} else { // for all other rows
				for (int col = 0; col < matrix.length; col++) { // iterate over columns of matrix
					if (col == columnRemoved) // when loop comes to desired column to be removed
						colShift = 1; // don't do anything with this column, and indicate column has been reached
					else // for all other columns
						minor[row - rowShift][col - colShift] = matrix[row][col];
					/*
					 * entry in minor set to entry in passed-in matrix, entries
					 * shifted up or left if removed row/column has been reached
					 */
				}
			}
		}
		return minor; // return matrix with row and column removed
	}

	public static double[] getRow(double[][] matrix, int rowIndex) { // gets row of matrix as vector
		return matrix[rowIndex]; // return row of matrix at desired row
	}

	public static double[] getColumn(double[][] matrix, int columnIndex) { // gets column of matrix as vector
		double[] columnArray = new double[matrix.length]; // init 1D array to hold column data
		for (int i = 0; i < matrix.length; i++) // iterate over rows
			columnArray[i] = matrix[i][columnIndex]; // fills column vector with entries of desired column of matrix
		return columnArray; // return column vector
	}

	public static void printMatrix(double[][] matrix) { // prints elements of 2D array 
		StringBuffer buff = new StringBuffer(); // inits new StringBuffer, allows appending entries
		for (int i = 0; i < matrix.length; i++) { // iterate over rows of passed matrix
			for (int j = 0; j < matrix[i].length; j++) { // iterate over columns of passed matrix
				buff.append(Math.round(matrix[i][j] * 10000.) / 10000.); // append entry to buff, rounded to 4 decimal places
				buff.append(" "); // add space between entries
			}
		}
		System.out.println(buff.toString()); // print matrix entries
	}
}
