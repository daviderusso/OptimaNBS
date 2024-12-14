package org.example;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Represents an instance of a spatial planning problem, including grid dimensions,
 * green infrastructure types, urban challenges, pre-existing and forbidden areas,
 * and scenario-specific data.
 */
public class Instance {
    // Dimensions of the grid (Width and Height).
    public int W;
    public int H;

    // List of green infrastructure types (e.g., Green Roof, Green Wall).
    public String[] GreenType;

    // List of urban challenges (e.g., Heat Island, Air Quality).
    public String[] UrbanChallenges;

    // Maps each GreenType to a list of points representing pre-existing green areas.
    public HashMap<String, ArrayList<Cells_Tile>> PreExistent;

    // Array storing the count of pre-existing areas for each GreenType.
    public int[] numberOfPreExistent;

    // Maps each GreenType to a boolean matrix indicating pre-existing green areas on the grid.
    public HashMap<String, boolean[][]> PreExistentMap;

    // Maps each GreenType to a list of points representing forbidden green areas.
    public HashMap<String, ArrayList<Cells_Tile>> Forbidden;

    // Array storing the count of forbidden areas for each GreenType.
    public int[] numberOfForbid;

    // Maps each GreenType to a boolean matrix indicating forbidden green areas on the grid.
    public HashMap<String, boolean[][]> ForbiddenMap;

    // 3D matrix storing Urban Challenges values for each urban challenge across the grid.
    public double[][][] A;

    // Maps each UrbanChallenge to its minimum and maximum values.
    public HashMap<String, double[]> BoundUCs;

    /**
     * Constructor for initializing an instance with specified attributes.
     * Calculates the minimum and maximum values for each urban challenge.
     *
     * @param w                Width of the grid.
     * @param h                Height of the grid.
     * @param greenType        Array of green infrastructure types.
     * @param urbanChallenges  Array of urban challenges.
     * @param preExistentMap   Boolean matrices indicating pre-existing green areas.
     * @param preExistent      Lists of points for pre-existing green areas.
     * @param countPreExistent Array of counts for pre-existing green areas.
     * @param forbiddenMap     Boolean matrices indicating forbidden green areas.
     * @param forbidden        Lists of points for forbidden green areas.
     * @param countForbid      Array of counts for forbidden green areas.
     * @param a                3D matrix of values for each urban challenge.
     */
    public Instance(int w, int h, String[] greenType, String[] urbanChallenges, HashMap<String, boolean[][]> preExistentMap,
                    HashMap<String, ArrayList<Cells_Tile>> preExistent, int[] countPreExistent, HashMap<String, boolean[][]> forbiddenMap,
                    HashMap<String, ArrayList<Cells_Tile>> forbidden, int[] countForbid, double[][][] a) {
        W = w;
        H = h;
        GreenType = greenType;
        PreExistentMap = preExistentMap;
        PreExistent = preExistent;
        numberOfPreExistent = countPreExistent;
        ForbiddenMap = forbiddenMap;
        Forbidden = forbidden;
        numberOfForbid = countForbid;
        A = a;
        UrbanChallenges = urbanChallenges;
        BoundUCs = new HashMap<>();

        // Calculate min and max values for each UrbanChallenge.
        for (int b = 0; b < urbanChallenges.length; b++) {
            double maxB = 0.0;
            double minB = 0.0;
            for (int i = 0; i < W; i++) {
                for (int j = 0; j < H; j++) {
                    if (A[b][i][j] < minB) {
                        minB = A[b][i][j];
                    }
                    if (A[b][i][j] > maxB) {
                        maxB = A[b][i][j];
                    }
                }
            }
            BoundUCs.put(urbanChallenges[b], new double[]{minB, maxB});
        }
    }

    /**
     * Default constructor for creating an empty instance with default values.
     */
    public Instance() {
        W = 1;
        H = 1;
        GreenType = new String[1];
        UrbanChallenges = new String[1];
        PreExistentMap = new HashMap<>();
        PreExistent = new HashMap<>();
        ForbiddenMap = new HashMap<>();
        Forbidden = new HashMap<>();
        A = new double[1][1][1];
    }

    /**
     * Checks if a specific cell (i, j) for a given green type is pre-existing.
     *
     * @param t The index of the GreenType.
     * @param i The row index in the grid.
     * @param j The column index in the grid.
     * @return True if the cell is pre-existing, false otherwise.
     */
    public boolean checkPreExist(int t, int i, int j) {
        return this.PreExistentMap.get(this.GreenType[t])[i][j];
    }

    /**
     * Prints detailed information about the instance, including grid dimensions,
     * green infrastructure types, urban challenges matrices.
     */
    public void printInstance() {
        System.out.println("---------------------------------------------");
        System.out.println(this.W + "x" + this.H);
        System.out.println("GreenType: " + Arrays.toString(this.GreenType));
        System.out.println("UrbanChallenges: " + Arrays.toString(this.UrbanChallenges));
        System.out.println();

        System.out.println("PreExistent Green areas: ");
        for (int g = 0; g < this.GreenType.length; g++) {
            System.out.println(this.GreenType[g] + ": " + this.numberOfPreExistent[g]);
            System.out.println(this.GreenType[g] + ": " + this.PreExistent.get(this.GreenType[g]));
        }
        System.out.println();

        System.out.println("Forbidden Green areas: ");
        for (int g = 0; g < this.GreenType.length; g++) {
            System.out.println(this.GreenType[g] + ": " + this.numberOfForbid[g]);
            System.out.println(this.GreenType[g] + ": " + this.Forbidden.get(this.GreenType[g]));
        }
        System.out.println();

        for (int a = 0; a < this.UrbanChallenges.length; a++) {
            System.out.println("Urban Challenges matrix for: " + this.UrbanChallenges[a]);
            for (int i = 0; i < this.W; i++) {
                for (int j = 0; j < this.H; j++) {
                    System.out.print(this.A[a][i][j] + "\t");
                }
                System.out.println();
            }
            System.out.println();
        }
        System.out.println("---------------------------------------------");
    }

    /**
     * Prints a summary of the instance, including green types and urban challenges.
     */
    public void printInstanceInfo() {
        System.out.println("---------------------------------------------");
        System.out.println(this.W + "x" + this.H);
        System.out.println("GreenType: " + Arrays.toString(this.GreenType));
        System.out.println("UrbanChallenges: " + Arrays.toString(this.UrbanChallenges));
        System.out.println("PreExistent Green areas: ");
        for (int g = 0; g < this.GreenType.length; g++) {
            System.out.println(this.GreenType[g] + ": " + this.numberOfPreExistent[g]);
        }
        System.out.println();

        System.out.println("Forbidden Green areas: ");
        for (int g = 0; g < this.GreenType.length; g++) {
            System.out.println(this.GreenType[g] + ": " + this.numberOfForbid[g]);
        }
        System.out.println("---------------------------------------------");
    }
}
