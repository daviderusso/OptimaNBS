package org.example;

/**
 * Represents a kernel structure associated with a specific GreenType and UrbanChallenge.
 * It includes the kernel's dimensions and its data matrix.
 */
public class Kernel {
    // The type of green infrastructure (e.g., Green Roof, Green Wall).
    public String GreenType;

    // The urban challenge this kernel addresses (e.g., Heat Island, Air Quality).
    public String UrbanChallenge;

    // A 2D array representing the kernel values.
    public double[][] k;

    // Width (number of columns) of the kernel matrix.
    public int w;

    // Height (number of rows) of the kernel matrix.
    public int h;

    /**
     * Constructs a Kernel object with the specified attributes.
     *
     * @param greenType      The type of green infrastructure.
     * @param urbanChallenge The urban challenge this kernel addresses.
     * @param w              The width of the kernel matrix.
     * @param h              The height of the kernel matrix.
     * @param k              A 2D array representing the kernel values.
     */
    public Kernel(String greenType, String urbanChallenge, int w, int h, double[][] k) {
        this.GreenType = greenType;
        this.UrbanChallenge = urbanChallenge;
        this.k = k;
        this.w = w;
        this.h = h;
    }

    /**
     * Prints the kernel details, including its metadata and matrix values.
     * The matrix is displayed row by row.
     */
    public void printKernel() {
        // Print header with kernel metadata.
        System.out.println("---------------------------------------------");
        System.out.println("UrbanChallenge: " + this.UrbanChallenge +
                " - GreenType: " + this.GreenType +
                " - Dimensions: " + this.w + "x" + this.h);
        System.out.println();

        // Print the kernel matrix.
        for (int i = 0; i < this.w; i++) {
            for (int j = 0; j < this.h; j++) {
                System.out.print(this.k[i][j] + "\t"); // Print each value in the row.
            }
            System.out.println(); // New line after each row.
        }
        System.out.println("---------------------------------------------");
    }
}
