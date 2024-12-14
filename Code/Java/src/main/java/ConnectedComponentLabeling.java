/**
 * This class implements the Connected Component Labeling (CCL) algorithm.
 * The algorithm identifies and labels distinct connected components in a binary matrix.
 */
public class ConnectedComponentLabeling {

    // Define the 8 possible neighbors (diagonal, horizontal, and vertical) for a cell.
    private static final int[][] NEIGHBORS = {
            {-1, -1}, {-1, 0}, {-1, 1},
            {0, -1}, {0, 1},
            {1, -1}, {1, 0}, {1, 1}
    };

    // Define the 4 possible neighbors (horizontal and vertical) for a cell.
    private static final int[][] NEIGHBORS_SMALL = {
            {-1, 0},  // Up
            {0, -1},  // Left
            {0, 1},   // Right
            {1, 0}    // Down
    };

    /**
     * Labels connected components in a binary matrix using depth-first search (DFS).
     * Each connected group of '1's in the binary matrix is assigned a unique label.
     *
     * @param binaryMatrix A 2D binary matrix (values are 0 or 1).
     * @return A 2D matrix with the same dimensions as the input, where each connected
     * component of '1's is labeled with a unique integer (starting from 1).
     */
    public static int[][] labelConnectedComponents(int[][] binaryMatrix) {
        int rows = binaryMatrix.length; // Number of rows in the matrix
        int cols = binaryMatrix[0].length; // Number of columns in the matrix
        int[][] labels = new int[rows][cols]; // Output matrix for labeled components
        int currentLabel = 0; // Current label to assign to a new connected component

        // Iterate over each cell in the matrix.
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                // If the cell is part of a component (value = 1) and has not been labeled.
                if (binaryMatrix[i][j] == 1 && labels[i][j] == 0) {
                    currentLabel++; // Increment the label for a new component.
                    // Perform DFS to label the entire connected component.
                    dfs(binaryMatrix, labels, i, j, currentLabel);
                }
            }
        }

        return labels; // Return the labeled matrix.
    }

    /**
     * Depth-first search (DFS) helper function to label all cells in a connected component.
     *
     * @param binaryMatrix The input binary matrix.
     * @param labels       The output matrix for storing component labels.
     * @param x            The row index of the current cell.
     * @param y            The column index of the current cell.
     * @param currentLabel The label to assign to the current connected component.
     */
    private static void dfs(int[][] binaryMatrix, int[][] labels, int x, int y, int currentLabel) {
        int rows = binaryMatrix.length; // Number of rows in the matrix.
        int cols = binaryMatrix[0].length; // Number of columns in the matrix.

        // Assign the current label to the cell.
        labels[x][y] = currentLabel;

        // Explore all neighbors of the current cell (using the NEIGHBORS_SMALL array).
        for (int[] neighbor : NEIGHBORS_SMALL) {
            int nx = x + neighbor[0]; // Row index of the neighbor.
            int ny = y + neighbor[1]; // Column index of the neighbor.

            // Check if the neighbor is within bounds, part of the component (value = 1),
            // and not already labeled.
            if (nx >= 0 && ny >= 0 && nx < rows && ny < cols &&
                    binaryMatrix[nx][ny] == 1 && labels[nx][ny] == 0) {
                // Recursively label the connected component.
                dfs(binaryMatrix, labels, nx, ny, currentLabel);
            }
        }
    }
}
