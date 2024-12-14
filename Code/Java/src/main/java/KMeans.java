import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * This class implements the K-Means clustering algorithm with K-Means++ initialization.
 * The algorithm partitions a set of points into K clusters, minimizing the variance within each cluster.
 */
public class KMeans {
    private int k; // Number of clusters.
    private static double tolerance = 1e-5; // Convergence tolerance for centroid movement.
    private List<PointKMeans> points; // List of points to be clustered.
    private List<Cluster> clusters; // List of clusters.

    /**
     * Constructor to initialize the KMeans instance.
     *
     * @param k      The number of clusters.
     * @param points The list of points to be clustered.
     */
    public KMeans(int k, List<PointKMeans> points) {
        this.k = k;
        this.points = points;
        this.clusters = new ArrayList<>();
    }

    /**
     * Performs the clustering using the K-Means algorithm.
     */
    public void cluster() {
        // Step 1: Initialize the clusters using K-Means++.
        initializeClustersKMeansPlusPlus();

        boolean converged = false;
        while (!converged) {
            // Step 2: Assign points to the nearest cluster.
            assignPoints();

            // Step 3: Update the centroids of the clusters.
            List<PointKMeans> prevCentroids = updateCentroids();

            // Step 4: Check if the centroids have converged (movement < tolerance).
            converged = true;
            for (int i = 0; i < k; i++) {
                if (clusters.get(i).centroid.distance(prevCentroids.get(i)) > tolerance) {
                    converged = false;
                    break;
                }
            }
        }
    }

    /**
     * Initializes the clusters using the K-Means++ initialization algorithm.
     */
    private void initializeClustersKMeansPlusPlus() {
        Random random = new Random();
        // Choose the first centroid randomly from the points.
        PointKMeans firstCentroid = points.get(random.nextInt(points.size()));
        clusters.add(new Cluster(firstCentroid));

        // Choose the remaining centroids based on distance-weighted probability.
        for (int i = 1; i < k; i++) {
            double[] distances = new double[points.size()];
            double sumDistanceSquared = 0;

            // Calculate the squared distance of each point to the nearest centroid.
            for (int j = 0; j < points.size(); j++) {
                double minDistanceSquared = Double.MAX_VALUE;
                for (Cluster cluster : clusters) {
                    double distanceSquared = points.get(j).distanceSquared(cluster.centroid);
                    minDistanceSquared = Math.min(minDistanceSquared, distanceSquared);
                }
                distances[j] = minDistanceSquared;
                sumDistanceSquared += minDistanceSquared;
            }

            // Choose the next centroid with probability proportional to squared distance.
            double randomValue = random.nextDouble() * sumDistanceSquared;
            double cumulativeProbability = 0;
            for (int j = 0; j < points.size(); j++) {
                cumulativeProbability += distances[j];
                if (cumulativeProbability >= randomValue) {
                    clusters.add(new Cluster(points.get(j)));
                    break;
                }
            }
        }
    }

    /**
     * Assigns each point to the nearest cluster based on the current centroids.
     */
    private void assignPoints() {
        // Clear all points from each cluster before reassignment.
        for (Cluster cluster : clusters) {
            cluster.points.clear();
        }

        // Assign each point to the nearest cluster.
        for (PointKMeans point : points) {
            Cluster closestCluster = null;
            double minDistance = Double.MAX_VALUE;

            for (Cluster cluster : clusters) {
                double distance = point.distance(cluster.centroid);
                if (distance < minDistance) {
                    minDistance = distance;
                    closestCluster = cluster;
                }
            }

            closestCluster.points.add(point); // Add the point to the closest cluster.
        }
    }

    /**
     * Updates the centroids of the clusters and returns the previous centroids.
     *
     * @return A list of the previous centroids.
     */
    private List<PointKMeans> updateCentroids() {
        List<PointKMeans> prevCentroids = new ArrayList<>();
        for (Cluster cluster : clusters) {
            PointKMeans centroid = cluster.calculateCentroid();
            prevCentroids.add(cluster.centroid); // Store the previous centroid.
            cluster.centroid = centroid; // Update to the new centroid.
        }
        return prevCentroids;
    }

    /**
     * Returns the list of clusters after clustering is complete.
     *
     * @return A list of clusters.
     */
    public List<Cluster> getClusters() {
        return clusters;
    }

    public static void main(String[] args) {
        // Example usage of KMeans.
        List<PointKMeans> points = new ArrayList<>();
        points.add(new PointKMeans(1, 2));
        points.add(new PointKMeans(1, 3));
        points.add(new PointKMeans(2, 2));
        points.add(new PointKMeans(5, 6));
        points.add(new PointKMeans(6, 5));
        points.add(new PointKMeans(7, 6));

        KMeans kmeans = new KMeans(2, points);
        kmeans.cluster();

        List<Cluster> clusters = kmeans.getClusters();
        for (int i = 0; i < clusters.size(); i++) {
            System.out.println("Cluster " + (i + 1) + ": " + clusters.get(i).points);
        }
    }
}

/**
 * Represents a point in 2D space.
 */
class PointKMeans {
    double x; // X-coordinate.
    double y; // Y-coordinate.

    public PointKMeans(double x, double y) {
        this.x = x;
        this.y = y;
    }

    /**
     * Calculates the Euclidean distance to another point.
     *
     * @param other The other point.
     * @return The Euclidean distance.
     */
    public double distance(PointKMeans other) {
        return Math.sqrt(distanceSquared(other));
    }

    /**
     * Calculates the squared Euclidean distance to another point.
     *
     * @param other The other point.
     * @return The squared Euclidean distance.
     */
    public double distanceSquared(PointKMeans other) {
        return Math.pow(x - other.x, 2) + Math.pow(y - other.y, 2);
    }

    @Override
    public String toString() {
        return "(" + x + ", " + y + ")";
    }
}

/**
 * Represents a cluster with a centroid and a list of points.
 */
class Cluster {
    PointKMeans centroid; // The centroid of the cluster.
    List<PointKMeans> points; // The points in the cluster.

    public Cluster(PointKMeans centroid) {
        this.centroid = centroid;
        this.points = new ArrayList<>();
    }

    /**
     * Calculates the new centroid based on the points in the cluster.
     *
     * @return The new centroid.
     */
    public PointKMeans calculateCentroid() {
        double sumX = 0, sumY = 0;
        for (PointKMeans point : points) {
            sumX += point.x;
            sumY += point.y;
        }
        double newX = sumX / points.size();
        double newY = sumY / points.size();
        return new PointKMeans(newX, newY);
    }
}
