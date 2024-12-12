package org.example;

public class Kernel {
    public String GreenType;
    public String Benefit;
    public double[][] k;

    public int w;
    public int h;

    public Kernel(String greenType, String benefit, int w, int h, double[][] k) {
        this.GreenType = greenType;
        this.Benefit = benefit;
        this.k = k;
        this.w = w;
        this.h = h;
    }

    public Kernel() {
        this.GreenType = "T1";
        this.Benefit = "B1";
        this.k = new double[3][3];
        this.w = 3;
        this.h = 3;
    }

    public void printKernel() {
        System.out.println("---------------------------------------------");
        System.out.println("Benefit: " + this.Benefit + " - " + "GreenType: " + this.GreenType + " - " + this.w + "x" + this.h);
        System.out.println();
        for (int i = 0; i < this.w; i++) {
            for (int j = 0; j < this.h; j++) {
                System.out.print(this.k[i][j] + "\t");
            }
            System.out.println();
        }
        System.out.println("---------------------------------------------");
    }
}
