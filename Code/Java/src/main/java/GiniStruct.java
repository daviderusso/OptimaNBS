import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;

/**
 * Class representing the data structure for calculating the Gini coefficient.
 * It stores the necessary data related to fairness (f) and population (pop) for each grid cell
 * and provides a method to compute the Gini index based on this data.
 */
class GiniData {
    // Attributes to store data for a particular cell (i, j)
    public double f;  // Fairness value for the cell
    public double pop;  // Population value for the cell
    public double f_perc;  // Fairness percentage of the total fairness
    public double pop_perc;  // Population percentage of the total population
    public double more_fair_pop_perc;  // Cumulative population percentage for more fair allocation
    public int i;  // Row index of the grid cell
    public int j;  // Column index of the grid cell

    /**
     * Constructor to initialize the GiniData object.
     *
     * @param i the row index of the grid cell
     * @param j the column index of the grid cell
     * @param f fairness value for the cell
     * @param pop population value for the cell
     */
    public GiniData(int i, int j, double f, double pop) {
        this.i = i;
        this.j = j;
        this.f = f;
        this.pop = pop;
    }

    @Override
    public String toString() {
        return i + "-" + j + " - f_perc: " + f_perc; // Return a string representation of the object
    }
}

/**
 * Class for computing the Gini index for fairness and population distribution.
 * It stores the fairness and population data and calculates the Gini coefficient,
 * which measures the inequality of distribution across a grid.
 */
public class GiniStruct {
    // List of GiniData objects representing the fairness and population for each cell
    ArrayList<GiniData> dataList;

    // Total fairness and population across all grid cells
    double f_tot;
    double pop_tot;

    // The calculated Gini coefficient
    double gini;

    public GiniStruct() {
        // Constructor to initialize the GiniStruct object
    }

    /**
     * Computes the Gini index based on fairness and population distribution.
     * The Gini index is a measure of inequality, where 0 represents perfect equality
     * and 1 represents maximum inequality.
     *
     * @param instance the Instance object containing grid data
     * @param resX the results of the installation (used for computing fairness)
     * @param kernel the kernel data used in fairness calculation
     * @param idFair the index of the fairness factor being analyzed
     * @return the Gini coefficient calculated based on the provided data
     */
    public double computeGini(Instance instance, int[][][] resX, HashMap<String, Kernel> kernel, int idFair) {
        // Initialize total fairness and population
        f_tot = 0.0;
        pop_tot = 0.0;
        dataList = new ArrayList<>();

        // Iterate over each grid cell (i, j) and calculate fairness and population data
        for (int i = 0; i < instance.W; i++) {
            for (int j = 0; j < instance.H; j++) {
                // Calculate fairness for the cell and add GiniData object to the list
                double fairness = InstanceUtils.computeKernelFairness(resX, i, j, kernel, instance.GreenType) * instance.A[idFair][i][j];
                double population = instance.A[idFair][i][j];
                dataList.add(new GiniData(i, j, fairness, population));

                // Accumulate the total fairness and population
                f_tot += fairness;
                pop_tot += population;
            }
        }

        // Calculate the percentage of fairness and population for each grid cell
        int progressiveID = 0;
        for (int i = 0; i < instance.W; i++) {
            for (int j = 0; j < instance.H; j++) {
                if (dataList.get(progressiveID).i == i && dataList.get(progressiveID).j == j) {
                    dataList.get(progressiveID).f_perc = dataList.get(progressiveID).f / f_tot;
                    dataList.get(progressiveID).pop_perc = dataList.get(progressiveID).pop / pop_tot;
                    progressiveID++;
                } else {
                    // This block should ideally not execute if data is consistent
                    System.out.println("ERROR - GINI COMPUTATION");
                }
            }
        }

        // Sort the data by the fairness percentage (f_perc) in ascending order
        dataList.sort(new Comparator<GiniData>() {
            public int compare(GiniData o1, GiniData o2) {
                if (o1.f_perc < o2.f_perc)
                    return -1;
                else if (o1.f_perc > o2.f_perc)
                    return 1;
                else
                    return 0;
            }
        });

        // Compute the cumulative population for more fair allocation
        double cumulativePop = 0.0;
        for (int i = 0, l = dataList.size(); i < l; i++) {
            cumulativePop += dataList.get(i).pop_perc;
            dataList.get(i).more_fair_pop_perc = 1.0 - cumulativePop; // More fair population percentage
        }

        // Calculate the Gini coefficient using the weighted sum of fairness values
        int nData = dataList.size();
        double coef = 2.0 / nData;
        double con = (nData + 1.0) / nData;
        double weighted_sum = 0.0;
        for (int i = 0; i < nData; i++) {
            weighted_sum += (i + 1) * dataList.get(i).f; // Weighted sum of fairness values
        }

        // Compute and return the Gini index
        gini = coef * (weighted_sum / f_tot) - con;
        return gini;
    }
}
