package org.example;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class InstanceUtils {

    /**
     * Reads an instance from a JSON file located at the given path.
     * Parses the instance properties, such as dimensions, GreenType, UrbanChallenges,
     * and various data maps, and creates an Instance object.
     *
     * @param path the path to the JSON file
     * @return the parsed Instance object
     */
    public static Instance readInstance(String path) {
        Instance inst = null;
        JSONParser jsonParser = new JSONParser();
        try (FileReader reader = new FileReader(path)) {
            // Parse the JSON file
            JSONObject json = (JSONObject) jsonParser.parse(reader);

            // Parse dimensions
            int W = Integer.parseInt(json.get(Config.W).toString());
            int H = Integer.parseInt(json.get(Config.H).toString());

            // Parse GreenType array
            JSONArray GreenTypeRaw = (JSONArray) json.get(Config.GT);
            String[] GreenType = new String[GreenTypeRaw.size()];
            for (int i = 0; i < GreenType.length; i++) {
                GreenType[i] = GreenTypeRaw.get(i).toString();
            }

            // Parse UrbanChallenges array
            JSONArray UCsRaw = (JSONArray) json.get(Config.UC);
            String[] UC = new String[UCsRaw.size()];
            for (int i = 0, l1 = UCsRaw.size(); i < l1; i++) {
                UC[i] = UCsRaw.get(i).toString();
            }

            // Parse PreExistent data
            JSONObject PreExistentRaw = (JSONObject) json.get(Config.PreExist);
            HashMap<String, boolean[][]> preExistentMap = new HashMap<>();
            HashMap<String, ArrayList<Cells_Tile>> preExistent = new HashMap<>();
            int[] countPreExistent = new int[GreenType.length];
            Arrays.fill(countPreExistent, 0);

            for (int i = 0; i < GreenType.length; i++) {
                boolean[][] tempPreExist = new boolean[W][H];
                for (int w = 0; w < W; w++) {
                    Arrays.fill(tempPreExist[w], false);
                }

                ArrayList<Cells_Tile> tempPoint = new ArrayList<>();
                JSONArray pointsRaw = (JSONArray) PreExistentRaw.get(GreenType[i]);

                for (int j = 0, l = pointsRaw.size(); j < l; j++) {
                    JSONArray currPoint = (JSONArray) pointsRaw.get(j);
                    tempPreExist[Integer.parseInt(currPoint.get(0).toString())]
                            [Integer.parseInt(currPoint.get(1).toString())] = true;
                    tempPoint.add(new Cells_Tile(Integer.parseInt(currPoint.get(0).toString()),
                            Integer.parseInt(currPoint.get(1).toString())));
                    countPreExistent[i]++;
                }
                preExistentMap.put(GreenType[i], tempPreExist);
                preExistent.put(GreenType[i], tempPoint);
            }

            // Parse Forbidden data
            JSONObject ForbiddenRaw = (JSONObject) json.get(Config.Forbid);
            HashMap<String, boolean[][]> forbiddenMap = new HashMap<>();
            HashMap<String, ArrayList<Cells_Tile>> forbidden = new HashMap<>();
            int[] countForbid = new int[GreenType.length];
            Arrays.fill(countForbid, 0);

            for (int i = 0; i < GreenType.length; i++) {
                boolean[][] tempForbid = new boolean[W][H];
                for (int w = 0; w < W; w++) {
                    Arrays.fill(tempForbid[w], false);
                }

                ArrayList<Cells_Tile> tempPoint = new ArrayList<>();
                JSONArray pointsRaw = (JSONArray) ForbiddenRaw.get(GreenType[i]);

                for (int j = 0, l = pointsRaw.size(); j < l; j++) {
                    JSONArray currPoint = (JSONArray) pointsRaw.get(j);
                    tempForbid[Integer.parseInt(currPoint.get(0).toString())]
                            [Integer.parseInt(currPoint.get(1).toString())] = true;
                    tempPoint.add(new Cells_Tile(Integer.parseInt(currPoint.get(0).toString()),
                            Integer.parseInt(currPoint.get(1).toString())));
                    countForbid[i]++;
                }
                forbiddenMap.put(GreenType[i], tempForbid);
                forbidden.put(GreenType[i], tempPoint);
            }

            // Parse A matrix
            JSONObject ARaw = (JSONObject) json.get(Config.A);
            double[][][] A = new double[UC.length][W][H];
            for (int b = 0; b < UC.length; b++) {
                JSONArray currA = (JSONArray) ARaw.get(UC[b]);
                int currAid = 0;
                for (int w = 0; w < W; w++) {
                    for (int h = 0; h < H; h++) {
                        A[b][w][h] = Double.parseDouble(currA.get(currAid++).toString());
                    }
                }
            }

            // Create Instance object
            inst = new Instance(W, H, GreenType, UC, preExistentMap, preExistent, countPreExistent,
                    forbiddenMap, forbidden, countForbid, A);
        } catch (IOException | ParseException e) {
            e.printStackTrace();
        }

        return inst;
    }

    /**
     * Computes the maximum fairness factor (FF) for each Urban Challenge without using Theta.
     * The computation considers forbidden cells and adjacency conditions.
     *
     * @param instance         the current instance
     * @param A                the benefit matrix
     * @param mat_of_forbidden the forbidden matrix
     * @param kernelSet        the kernel set
     * @return a map of Urban Challenges to their maximum FF
     */
    public static HashMap<String, Double> computeMaxFFNoTheta(Instance instance, double[][][] A, int[][] mat_of_forbidden, HashMap<String, HashMap<String, Kernel>> kernelSet) {
        HashMap<String, Double> singleFF = new HashMap<>();
        double currentFF;

        for (int b = 0, l2 = instance.UrbanChallenges.length; b < l2; b++) {
            currentFF = 0.0;

            for (int i = 0; i < instance.W; i++) {
                for (int j = 0; j < instance.H; j++) {
                    if (currentFF < A[b][i][j]) {
                        if ((mat_of_forbidden[i][j] < instance.GreenType.length) ||
                                (checkAdjTile(instance, mat_of_forbidden, i, j, kernelSet))) {
                            currentFF = A[b][i][j];
                        }
                    }
                }
            }
            singleFF.put(instance.UrbanChallenges[b], currentFF);
        }

        return singleFF;
    }

    /**
     * Checks adjacency conditions for a given tile.
     * Determines if a tile is adjacent to a forbidden cell based on kernel conditions.
     *
     * @param instance         the current instance
     * @param mat_of_forbidden the forbidden matrix
     * @param i                the row index of the tile
     * @param j                the column index of the tile
     * @param kernelSet        the kernel set
     * @return true if the tile is adjacent to a forbidden cell, false otherwise
     */
    public static boolean checkAdjTile(Instance instance, int[][] mat_of_forbidden, int i, int j, HashMap<String, HashMap<String, Kernel>> kernelSet) {
        int W = mat_of_forbidden.length;
        int H = mat_of_forbidden[0].length;

        for (int uc = 0, l2 = instance.UrbanChallenges.length; uc < l2; uc++) {
            for (int t = 0; t < instance.GreenType.length; t++) {
                int sizeKw = kernelSet.get(instance.UrbanChallenges[uc]).get(instance.GreenType[t]).w;
                int sizeKh = kernelSet.get(instance.UrbanChallenges[uc]).get(instance.GreenType[t]).h;
                int a = sizeKw / 2;
                int b = sizeKh / 2;

                for (int ii = 0; ii < a + 1; ii++) {
                    if (i - (a - ii) >= 0) {
                        for (int jj = 0; jj < b + 1; jj++) {
                            if (j - (b - jj) >= 0 && mat_of_forbidden[i - (a - ii)][j - (b - jj)] < instance.GreenType.length) {
                                return true;
                            }
                        }
                        for (int jj = b + 1; jj < sizeKh; jj++) {
                            if (j + (jj - b) < H && mat_of_forbidden[i - (a - ii)][j + (jj - b)] < instance.GreenType.length) {
                                return true;
                            }
                        }
                    }
                }
                for (int ii = a + 1; ii < sizeKw; ii++) {
                    if (i + (ii - a) < W) {
                        for (int jj = 0; jj < b + 1; jj++) {
                            if (j - (b - jj) >= 0 && mat_of_forbidden[i + (ii - a)][j - (b - jj)] < instance.GreenType.length) {
                                return true;
                            }
                        }
                        for (int jj = b + 1; jj < sizeKh; jj++) {
                            if (j + (jj - b) < H && mat_of_forbidden[i + (ii - a)][j + (jj - b)] < instance.GreenType.length) {
                                return true;
                            }
                        }
                    }
                }
            }
        }
        return false;
    }

    /**
     * Computes the size of each cluster in a clustered map.
     *
     * @param clusteredMapPark the clustered map containing cluster IDs
     * @return a list of cluster sizes indexed by cluster ID
     */
    public static ArrayList<Integer> computeClusterSize(int[][] clusteredMapPark) {
        ArrayList<Integer> sizeCluster = new ArrayList<>();
        int W = clusteredMapPark.length;
        int H = clusteredMapPark[0].length;

        for (int w = 0; w < W; w++) {
            for (int h = 0; h < H; h++) {
                if (clusteredMapPark[w][h] != 0) {
                    if (sizeCluster.size() <= clusteredMapPark[w][h]) {
                        while (sizeCluster.size() <= clusteredMapPark[w][h]) {
                            sizeCluster.add(0);
                        }
                    }
                    sizeCluster.set(clusteredMapPark[w][h], sizeCluster.get(clusteredMapPark[w][h]) + 1);
                }
            }
        }
        return sizeCluster;
    }

    /**
     * Computes the fitness function (FF) for an urban challenge instance based on a given configuration.
     *
     * @param instance         the urban challenge instance containing data related to the challenge
     * @param config           the configuration object containing the challenge settings and parameters
     * @param A                the matrix containing challenge data for fitness computation
     * @param mat_of_forbidden a matrix indicating forbidden locations on the map
     * @param kernelSet        a map containing kernels for additional checks during fitness computation
     * @return a HashMap containing the fitness function values for each urban challenge key
     */
    public static HashMap<String, Double> computeOF(Instance instance, Config config, double[][][] A, int[][] mat_of_forbidden, HashMap<String, HashMap<String, Kernel>> kernelSet) {
        HashMap<String, Double> singleFF = new HashMap<>();
        double currentFF = 0.0;

        // Iterate over all keys in the UCMap of the config
        for (String key : config.UCMap.keySet()) {

            // Find the position of the current urban challenge key in the A matrix
            int currB = -1;
            for (int b = 0; b < instance.UrbanChallenges.length; b++) {
                if (Objects.equals(key, instance.UrbanChallenges[b])) {
                    currB = b;
                    break;
                }
            }

            // For each value associated with the key in the UCMap
            for (int bb = 0, l2 = config.UCMap.get(key).length; bb < l2; bb++) {

                // If the current element is 'Avg', compute the average fitness function
                if (Objects.equals(config.UCMap.get(key)[bb], Config.Avg)) {
                    currentFF = 0.0;
                    // Iterate over the grid cells (W x H)
                    for (int i = 0; i < instance.W; i++) {
                        for (int j = 0; j < instance.H; j++) {
                            // If the cell is not forbidden or satisfies adjacency conditions
                            if ((mat_of_forbidden[i][j] < instance.GreenType.length) || (checkAdjTile(instance, mat_of_forbidden, i, j, kernelSet))) {
                                currentFF += A[currB][i][j]; // Add the value from A to the sum
                            }
                        }
                    }
                    // Calculate the average and store the result
                    currentFF = currentFF / (instance.H * instance.W);
                    singleFF.put(key + "_" + config.UCMap.get(key)[bb], currentFF);
                }
                // If the current element is 'Max', compute the maximum fitness function
                else if (Objects.equals(config.UCMap.get(key)[bb], Config.Max)) {
                    currentFF = 0.0;
                    // Iterate over the grid cells (W x H)
                    for (int i = 0; i < instance.W; i++) {
                        for (int j = 0; j < instance.H; j++) {
                            if (currentFF < A[currB][i][j]) {
                                // If the value at A is higher and the cell is valid
                                if ((mat_of_forbidden[i][j] < instance.GreenType.length) || (checkAdjTile(instance, mat_of_forbidden, i, j, kernelSet))) {
                                    currentFF = A[currB][i][j]; // Update the maximum fitness function
                                }
                            }
                        }
                    }
                    // Store the result
                    singleFF.put(key + "_" + config.UCMap.get(key)[bb], currentFF);
                }
            }
        }
        return singleFF;
    }

    /**
     * Computes the Gini coefficient for an urban challenge instance based on a given result matrix.
     *
     * @param instance the urban challenge instance containing the challenge data
     * @param resX     the result matrix used to compute the Gini coefficient
     * @param kernel   a map containing kernel functions for additional checks during the computation
     * @return the computed Gini coefficient, a measure of inequality
     */
    public static double computeGiniCoefficient(Instance instance, int[][][] resX, HashMap<String, Kernel> kernel) {
        double giniCoeff = 0.0;
        int idFair = instance.UrbanChallenges.length - 1;

        // Compute the Gini coefficient using a helper structure
        GiniStruct struct = new GiniStruct();
        giniCoeff = struct.computeGini(instance, resX, kernel, idFair);

        // Ensure the Gini coefficient is within the valid range [0, 1]
        if (giniCoeff < 0.0) {
            giniCoeff = 0.0;
        } else if (giniCoeff > 1.0 || Double.isNaN(giniCoeff)) {
            giniCoeff = 1.0;
        }

        return giniCoeff;
    }

    /**
     * Computes the kernel fairness for a specific position in a grid based on the result matrix and kernel values.
     *
     * @param resX      the result matrix used for computing fairness
     * @param i         the row index in the grid
     * @param j         the column index in the grid
     * @param kernel    a map containing the kernels used for fairness computation
     * @param GreenType an array of green types, each representing a specific challenge type
     * @return the computed kernel fairness value for the specified grid position
     */
    public static double computeKernelFairness(int[][][] resX, int i, int j, HashMap<String, Kernel> kernel, String[] GreenType) {
        double val = 0.0; // Initialize the fairness value to 0.0
        double[][] K; // Variable to store the kernel matrix
        int a, b; // Variables to store the dimensions of the kernel
        int W = resX[0].length; // Width of the result matrix
        int H = resX[0][0].length; // Height of the result matrix

        // Iterate over all green types
        for (int t = 0; t < GreenType.length; t++) {
            K = kernel.get(GreenType[t]).k; // Get the kernel for the current green type
            a = K.length / 2; // Calculate the half-width of the kernel
            b = K[0].length / 2; // Calculate the half-height of the kernel

            // Apply the kernel to the result matrix for positions above and to the left of (i, j)
            for (int ii = 0; ii < a + 1; ii++) {
                if (i - (a - ii) >= 0) { // Check if the row index is valid
                    for (int jj = 0; jj < b + 1; jj++) {
                        if (j - (b - jj) >= 0) { // Check if the column index is valid
                            val += (K[ii][jj] * resX[t][i - (a - ii)][j - (b - jj)]); // Compute fairness contribution
                        }
                    }

                    // Apply the kernel to positions above and to the right of (i, j)
                    for (int jj = b + 1; jj < K.length; jj++) {
                        if (j + (jj - b) < H) {
                            val += (K[ii][jj] * resX[t][i - (a - ii)][j + (jj - b)]);
                        }
                    }
                }
            }

            // Apply the kernel to positions below and to the left of (i, j)
            for (int ii = a + 1; ii < K.length; ii++) {
                if (i + (ii - a) < W) {
                    for (int jj = 0; jj < b + 1; jj++) {
                        if (j - (b - jj) >= 0) {
                            val += (K[ii][jj] * resX[t][i + (ii - a)][j - (b - jj)]);
                        }
                    }

                    // Apply the kernel to positions below and to the right of (i, j)
                    for (int jj = b + 1; jj < K[0].length; jj++) {
                        if (j + (jj - b) < H) {
                            val += (K[ii][jj] * resX[t][i + (ii - a)][j + (jj - b)]);
                        }
                    }
                }
            }
        }
        return val; // Return the computed kernel fairness value
    }

    /**
     * Computes both the minimum and maximum fairness values for an urban challenge instance based on pre-existing data and kernel values.
     *
     * @param instance the urban challenge instance containing all relevant data, including GreenType, Forbidden areas, and kernel configurations
     * @param kernel   a map containing the kernels used for fairness computation
     * @return an array containing the minimum and maximum fairness values, where index 0 holds the minimum fairness and index 1 holds the maximum fairness
     */
    public static double[] computeBoundFairness(Instance instance, HashMap<String, Kernel> kernel) {
        double MinFairness = 0.0; // Initialize the minimum fairness value
        int idFair = instance.UrbanChallenges.length - 1; // Get the index for the fairness challenge

        // Build the pre-existing matrix (preExistentX) based on the Forbidden areas
        int[][][] preExistentX = new int[instance.GreenType.length][instance.W][instance.H];
        for (int i = 0; i < instance.GreenType.length; i++) {
            for (int j = 0, l = instance.PreExistent.get(instance.GreenType[i]).size(); j < l; j++) {
                preExistentX[i][instance.PreExistent.get(instance.GreenType[i]).get(j).r][instance.PreExistent.get(instance.GreenType[i]).get(j).c] = 1;
            }
        }

        // Compute the minimum fairness using pre-existing data
        for (int i = 0; i < instance.W; i++) {
            for (int j = 0; j < instance.H; j++) {
                MinFairness += computeKernelFairness(preExistentX, i, j, kernel, instance.GreenType) * instance.A[idFair][i][j];
            }
        }

        // Compute the maximum fairness
        double MaxFairness = 0.0;
        double maxVal = 0;
        int Gtmax = 0;
        for (int i = 0; i < instance.GreenType.length; i++) {
            int i1 = kernel.get(instance.GreenType[i]).w / 2; // Get half the width of the kernel
            int i2 = kernel.get(instance.GreenType[i]).h / 2; // Get half the height of the kernel
            double val = kernel.get(instance.GreenType[i]).k[i1][i2]; // Get the central kernel value
            if (val > maxVal) { // Find the kernel with the maximum value
                maxVal = val;
                Gtmax = i;
            }
        }

        // Build the forbidden map matrix (mat_of_forbidden) based on the Forbidden areas
        int[][] mat_of_forbidden = new int[instance.W][instance.H];
        int curr_i, curr_j;
        for (int t = 0, len = instance.GreenType.length; t < len; t++) {
            for (int f = 0, lenF = instance.Forbidden.get(instance.GreenType[t]).size(); f < lenF; f++) {
                curr_i = instance.Forbidden.get(instance.GreenType[t]).get(f).r;
                curr_j = instance.Forbidden.get(instance.GreenType[t]).get(f).c;
                mat_of_forbidden[curr_i][curr_j]++;
            }
        }

        // Build the ideal matrix (idealX) based on the forbidden map
        int[][][] idealX = new int[instance.GreenType.length][instance.W][instance.H];
        for (int i = 0; i < instance.W; i++) {
            for (int j = 0; j < instance.H; j++) {
                if (mat_of_forbidden[i][j] < instance.GreenType.length) {
                    idealX[Gtmax][i][j] = 1; // Mark the ideal positions for the green type with the maximum kernel value
                    for (int k = 0; k < instance.GreenType.length; k++) {
                        if (k != Gtmax) {
                            idealX[k][i][j] = 0; // Set other green types to 0
                        }
                    }
                } else {
                    idealX[Gtmax][i][j] = 0; // Set non-ideal positions to 0
                }
            }
        }

        // Compute the maximum fairness using the ideal matrix
        for (int i = 0; i < instance.W; i++) {
            for (int j = 0; j < instance.H; j++) {
                MaxFairness += computeKernelFairness(idealX, i, j, kernel, instance.GreenType) * instance.A[idFair][i][j];
            }
        }

        return new double[]{MinFairness, MaxFairness}; // Return both the minimum and maximum fairness values
    }

    /**
     * Computes a clustered map with connected component labeling (CCL) and applies size-based filtering and splitting for clusters.
     * The method labels connected components in a grid, filters clusters based on size constraints, and splits large clusters into smaller ones.
     *
     * @param instance the urban challenge instance containing the map and forbidden areas
     * @param config   the configuration object containing the minimum, maximum, and limit sizes for cluster filtering
     * @return a clustered map with labeled clusters that satisfy the specified size constraints
     */
    public static int[][] getParkClusterCCLWithBounds(Instance instance, Config config) {
        int[][] clusteredMap = new int[instance.W][instance.H]; // Initialize the clustered map with dimensions W x H
        int[][] tempClusteredMap = new int[instance.W][instance.H]; // Temporary map for storing cluster updates

        // Initialize all map cells as belonging to cluster 1
        for (int w = 0; w < instance.W; w++) {
            Arrays.fill(clusteredMap[w], 1);
        }

        try {
            // Set forbidden areas (defined in the configuration) to cluster 0
            for (int i = 0, l = instance.Forbidden.get(Config.ClusteredGT).size(); i < l; i++) {
                clusteredMap[instance.Forbidden.get(Config.ClusteredGT).get(i).r][instance.Forbidden.get(Config.ClusteredGT).get(i).c] = 0;
            }

            // Perform connected component labeling to group connected cells into clusters
            clusteredMap = ConnectedComponentLabeling.labelConnectedComponents(clusteredMap);

            // Compute the size of each cluster
            ArrayList<Integer> sizeCluster = InstanceUtils.computeClusterSize(clusteredMap);

            // Filter clusters based on size constraints: small clusters are removed, good clusters are preserved, and large clusters are split
            for (int i = 1, l = sizeCluster.size(); i < l; i++) {
                if (sizeCluster.get(i) < Config.minClusterSize) {
                    // Remove small clusters by setting their cells to 0 in the temporary map
                    for (int w = 0; w < instance.W; w++) {
                        for (int h = 0; h < instance.H; h++) {
                            if (clusteredMap[w][h] == i) {
                                tempClusteredMap[w][h] = 0;
                            }
                        }
                    }
                } else if (sizeCluster.get(i) < Config.maxClusterSize) {
                    // Preserve good-sized clusters in the temporary map
                    for (int w = 0; w < instance.W; w++) {
                        for (int h = 0; h < instance.H; h++) {
                            if (clusteredMap[w][h] == i) {
                                tempClusteredMap[w][h] = clusteredMap[w][h];
                            }
                        }
                    }
                } else if (sizeCluster.get(i) > Config.maxLimitSize) {
                    // Split large clusters by applying KMeans to remove branching tiles
                    List<PointKMeans> points = new ArrayList<>();
                    for (int w = 0; w < instance.W; w++) {
                        for (int h = 0; h < instance.H; h++) {
                            if (clusteredMap[w][h] == i) {
                                points.add(new PointKMeans(w, h));
                            }
                        }
                    }
                    int nCluster = points.size() / Config.maxClusterSize; // Determine the number of clusters for splitting
                    KMeans kmeans = new KMeans(nCluster, points); // Apply KMeans clustering
                    kmeans.cluster();
                    List<Cluster> clusters = kmeans.getClusters();
                    int newClusterId = sizeCluster.size();

                    // Assign new cluster IDs after splitting
                    for (int j = 0, l1 = clusters.size(); j < l1; j++) {
                        for (int k = 0, l2 = clusters.get(j).points.size(); k < l2; k++) {
                            tempClusteredMap[(int) clusters.get(j).points.get(k).x][(int) clusters.get(j).points.get(k).y] = newClusterId;
                        }
                        newClusterId++;
                    }
                }
            }

            // Iteratively split clusters that exceed the size limit
            boolean overSize;
            int test = 0;
            do {
                overSize = false;
                sizeCluster = InstanceUtils.computeClusterSize(tempClusteredMap); // Recompute cluster sizes
                for (int i = 1, l = sizeCluster.size(); i < l; i++) {
                    if (sizeCluster.get(i) > Config.maxLimitSize) {
                        overSize = true;
                        // Split the oversized clusters again using KMeans
                        List<PointKMeans> points = new ArrayList<>();
                        for (int w = 0; w < instance.W; w++) {
                            for (int h = 0; h < instance.H; h++) {
                                if (tempClusteredMap[w][h] == i) {
                                    points.add(new PointKMeans(w, h));
                                }
                            }
                        }
                        int nCluster = points.size() / Config.maxClusterSize;
                        KMeans kmeans = new KMeans(nCluster, points);
                        kmeans.cluster();
                        List<Cluster> clusters = kmeans.getClusters();
                        int newClusterId = sizeCluster.size();

                        // Update the cluster IDs after splitting
                        for (int j = 0, l1 = clusters.size(); j < l1; j++) {
                            for (int k = 0, l2 = clusters.get(j).points.size(); k < l2; k++) {
                                tempClusteredMap[(int) clusters.get(j).points.get(k).x][(int) clusters.get(j).points.get(k).y] = newClusterId;
                            }
                            newClusterId++;
                        }
                    }
                }
                test++;
            } while (overSize && test < 50); // Limit to 50 iterations to avoid infinite loop

            // Normalize the clusters by removing empty clusters
            sizeCluster = InstanceUtils.computeClusterSize(tempClusteredMap);
            int currID = 1;
            for (int i = 1, l = sizeCluster.size(); i < l; i++) {
                if (sizeCluster.get(i) != 0) {
                    // Assign new cluster IDs for non-empty clusters
                    for (int w = 0; w < instance.W; w++) {
                        for (int h = 0; h < instance.H; h++) {
                            if (tempClusteredMap[w][h] == i) {
                                clusteredMap[w][h] = currID;
                            }
                        }
                    }
                    currID++;
                }
            }

            sizeCluster = InstanceUtils.computeClusterSize(clusteredMap); // Recompute the cluster sizes after normalization

        } catch (Exception e) {
            System.out.println("CLUSTERING ERROR");
        }
        return clusteredMap; // Return the final clustered map
    }

    /**
     * Normalizes the data in the given Instance object.
     * This method normalizes the attributes (A) of the instance based on predefined bounds,
     * and also normalizes the cost and budget based on the configuration.
     *
     * @param instance the Instance object to be normalized
     * @param config   the configuration object containing the cost and budget details
     * @return the normalized instance with updated attributes and costs
     */
    public static Instance normalizeData(Instance instance, Config config) {
        // Normalize the attributes A of the instance using the bounds for each urban challenge
        for (int b = 0; b < instance.UrbanChallenges.length; b++) {
            for (int w = 0; w < instance.W; w++) {
                for (int h = 0; h < instance.H; h++) {
                    // Normalize the value of A[b][w][h] to be between 0 and 1 using the bounds
                    instance.A[b][w][h] = (instance.A[b][w][h] - instance.BoundUCs.get(instance.UrbanChallenges[b])[0])
                            / (instance.BoundUCs.get(instance.UrbanChallenges[b])[1] - instance.BoundUCs.get(instance.UrbanChallenges[b])[0]);
                    // Ensure values remain within the [0, 1] range, correcting any out-of-bound values
                    if (instance.A[b][w][h] > 1.0 || instance.A[b][w][h] < 0.0 || Double.isNaN(instance.A[b][w][h])) {
                        if (instance.A[b][w][h] > 1.0) {
                            instance.A[b][w][h] = 1.0;
                        } else {
                            instance.A[b][w][h] = 0.0;
                        }
                    }
                }
            }
        }

        // Normalize the cost for each GreenType based on the maximum cost
        double c_max = 0.0;
        double c_min = 0.0;
        for (int t = 0; t < instance.GreenType.length; t++) {
            if (c_max < config.cost.get(instance.GreenType[t])) {
                c_max = config.cost.get(instance.GreenType[t]);
            }
        }
        c_max = c_max * instance.H * instance.W;

        // Store the normalized cost in the configuration object
        for (int t = 0; t < instance.GreenType.length; t++) {
            config.normalizedCost.put(instance.GreenType[t], ((config.cost.get(instance.GreenType[t]) - c_min) / (c_max - c_min)));
        }

        // Normalize the budget based on the min and max cost
        config.NormalizedBudget = (config.Budget - c_min) / (c_max - c_min);
        if (config.NormalizedBudget > 1.0) {
            config.NormalizedBudget = 1.0;
        } else if (config.NormalizedBudget < 0.0) {
            config.NormalizedBudget = 0.0;
        }

        return instance; // Return the normalized instance
    }

    /**
     * Denormalizes the data in the given Instance object by reversing the normalization process.
     * The method restores the original values for the attributes (A) using the stored bounds.
     *
     * @param instance the Instance object to be denormalized
     * @return the denormalized instance with original attribute values
     */
    public static Instance denormalizeData(Instance instance) {
        // Restore original values for each attribute A by reversing the normalization
        for (int b = 0; b < instance.UrbanChallenges.length; b++) {
            for (int w = 0; w < instance.W; w++) {
                for (int h = 0; h < instance.H; h++) {
                    instance.A[b][w][h] = (instance.A[b][w][h] * ((instance.BoundUCs.get(instance.UrbanChallenges[b])[1] - instance.BoundUCs.get(instance.UrbanChallenges[b])[0])))
                            + instance.BoundUCs.get(instance.UrbanChallenges[b])[0];
                }
            }
        }
        return instance; // Return the denormalized instance
    }

    /**
     * Denormalizes the given 3D array of attribute values A by reversing the normalization process.
     * This method is used when you want to apply denormalization to a copy of the attribute values without modifying the original data.
     *
     * @param instance the Instance object containing the bounds for denormalization
     * @param A        the 3D array of normalized attribute values to be denormalized
     * @return the denormalized 3D array of attribute values
     */
    public static double[][][] denormalizeData(Instance instance, double[][][] A) {
        // Restore original values for each attribute A in the 3D array by reversing the normalization
        for (int b = 0; b < instance.UrbanChallenges.length; b++) {
            for (int w = 0; w < instance.W; w++) {
                for (int h = 0; h < instance.H; h++) {
                    A[b][w][h] = (A[b][w][h] * ((instance.BoundUCs.get(instance.UrbanChallenges[b])[1] - instance.BoundUCs.get(instance.UrbanChallenges[b])[0])))
                            + instance.BoundUCs.get(instance.UrbanChallenges[b])[0];
                }
            }
        }
        return A; // Return the denormalized 3D array
    }

    /**
     * Normalizes the given 3D array of attribute values A to the range [0, 1] using the bounds stored in the instance.
     * This method is used when you want to normalize a separate copy of the attribute values without modifying the original data.
     *
     * @param instance the Instance object containing the bounds for normalization
     * @param A        the 3D array of attribute values to be normalized
     * @return the normalized 3D array of attribute values
     */
    public static double[][][] normalizeData(Instance instance, double[][][] A) {
        // Normalize each attribute in the 3D array A using the bounds for each urban challenge
        for (int b = 0; b < instance.UrbanChallenges.length; b++) {
            for (int w = 0; w < instance.W; w++) {
                for (int h = 0; h < instance.H; h++) {
                    A[b][w][h] = (A[b][w][h] - instance.BoundUCs.get(instance.UrbanChallenges[b])[0])
                            / (instance.BoundUCs.get(instance.UrbanChallenges[b])[1] - instance.BoundUCs.get(instance.UrbanChallenges[b])[0]);
                    // Check if the normalization went out of bounds (should be [0, 1])
                    if (A[b][w][h] > 1) {
                        System.out.println("ERROR: normalization data failed");
                    }
                }
            }
        }
        return A; // Return the normalized 3D array
    }


}