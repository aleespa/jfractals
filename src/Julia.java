import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import static java.lang.Math.log;
import static java.lang.Math.sqrt;

public class Julia {

    static final int M = 4000;
    static final int N = 4000;
    static final double S = 2500.0;
    static final int MAX_ITER = 256;

    public static void main(String[] args) {
        final double real = Double.parseDouble(args[0]);
        final double imag = Double.parseDouble(args[1]);
        final String filename = args[2];
        double[][] iterations = new double[N][M];

        // Compute Mandelbrot set
        juliaSet(iterations, real, imag);

        // Save results to a file
        saveToFile(iterations, filename);

        // Print a sample value for verification
        System.out.println("Iterations at center: " + iterations[N / 2][M / 2]);
    }

    public static void juliaSet(double[][] iterations,
                                double real,
                                double imag) {
        double[] xGrid = linspace(-M / S, M / S, M);
        double[] yGrid = linspace(-N / S, N / S, N);

        java.util.stream.IntStream.range(0, N).parallel().forEach(i -> {
            for (int j = 0; j < M; j++) {
                double x = xGrid[j];
                double y = yGrid[i];
                double zx = x;
                double zy = y;
                double iter = 0;
                double lastMagnitude =0;

                for (int k = 0; k < MAX_ITER; k++) {
                    double zx2 = zx * zx;
                    double zy2 = zy * zy;

                    if (zx2 + zy2 > 1000) {
                        break;
                    }

                    double temp = zx2 - zy2 + real;
                    zy = 2.0 * zx * zy + imag;
                    zx = temp;

                    iter ++;
                    lastMagnitude = sqrt(zx2 + zy2);
                }

                if (lastMagnitude < 1){
                    iterations[i][j] = MAX_ITER;
                } else{
                    iterations[i][j] = round(iter + 1 - log(log(lastMagnitude)) / log(2),2);
                }
            }
        });
    }

    // Function to create a linearly spaced array
    public static double[] linspace(double start, double end, int num) {
        double[] result = new double[num];
        double step = (end - start) / (num - 1);
        for (int i = 0; i < num; i++) {
            result[i] = start + i * step;
        }
        return result;
    }

    public static void saveToFile(double[][] iterations, String filename) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filename))) {
            for (int i = 0; i < iterations.length; i++) {
                for (int j = 0; j < iterations[i].length; j++) {
                    writer.write(iterations[i][j] + (j == iterations[i].length - 1 ? "" : ","));
                }
                writer.newLine();
            }
            System.out.println("Mandelbrot set saved to " + filename);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static double round(double a, int k){
        return Math.round(a * Math.pow(10, k)) /  Math.pow(10, k);
    }
}
