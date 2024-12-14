package org.example;

/**
 * This class represents a cell in a tile-based grid.
 * Each cell is defined by its row (`r`) and column (`c`) indices.
 */
public class Cells_Tile {
    int r; // Row index of the cell.
    int c; // Column index of the cell.

    /**
     * Constructor to initialize a cell with specified row and column indices.
     *
     * @param r The row index of the cell.
     * @param c The column index of the cell.
     */
    public Cells_Tile(int r, int c) {
        this.r = r;
        this.c = c;
    }

    /**
     * Overrides the `toString` method to provide a string representation of the cell.
     * The format is "(row, column)".
     *
     * @return A string representation of the cell.
     */
    @Override
    public String toString() {
        return "(" + r + " , " + c + ")";
    }
}
