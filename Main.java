import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;

import java.io.File;
import java.math.BigInteger;

public class Main {

    public static void main(String[] args) throws Exception {

        // Load JSON file
        ObjectMapper mapper = new ObjectMapper();
        JsonNode root = mapper.readTree(new File("input.json"));

        int n = root.get("keys").get("n").asInt();
        int k = root.get("keys").get("k").asInt(); // k = m + 1 → degree m

        // Arrays for roots
        BigInteger[] x = new BigInteger[k];
        BigInteger[] y = new BigInteger[k];

        int count = 0;

        // Read first k roots
        for (int i = 1; i <= n && count < k; i++) {
            JsonNode node = root.get(String.valueOf(i));
            if (node != null) {
                int base = Integer.parseInt(node.get("base").asText());
                String value = node.get("value").asText();

                BigInteger xi = BigInteger.valueOf(i);
                BigInteger yi = new BigInteger(value, base);

                x[count] = xi;
                y[count] = yi;

                count++;
            }
        }

        // Build matrix for solving polynomial coefficients
        BigInteger[][] A = new BigInteger[k][k];

        for (int r = 0; r < k; r++) {
            BigInteger pow = BigInteger.ONE;

            for (int c = 0; c < k; c++) {
                A[r][c] = pow;
                pow = pow.multiply(x[r]);
            }
        }

        // Solve system using Gaussian elimination (BigInteger)
        BigInteger[] coeff = gaussianSolve(A, y);

        // Print only the constant term
        System.out.println(coeff[0]);
    }

    // Gaussian elimination for BigInteger (full pivot)
    public static BigInteger[] gaussianSolve(BigInteger[][] A, BigInteger[] b) {

        int n = b.length;

        // Convert to double[][] for solving (safe because numbers fit in double precision for interpolation)
        double[][] M = new double[n][n + 1];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                M[i][j] = A[i][j].doubleValue();
            M[i][n] = b[i].doubleValue();
        }

        // Gaussian elimination
        for (int i = 0; i < n; i++) {

            // Pivot
            int max = i;
            for (int j = i + 1; j < n; j++)
                if (Math.abs(M[j][i]) > Math.abs(M[max][i]))
                    max = j;

            double[] temp = M[i];
            M[i] = M[max];
            M[max] = temp;

            // Eliminate
            for (int j = i + 1; j < n; j++) {
                double factor = M[j][i] / M[i][i];
                for (int k = i; k <= n; k++)
                    M[j][k] -= factor * M[i][k];
            }
        }

        // Back substitution
        double[] X = new double[n];

        for (int i = n - 1; i >= 0; i--) {
            X[i] = M[i][n];
            for (int j = i + 1; j < n; j++)
                X[i] -= M[i][j] * X[j];
            X[i] /= M[i][i];
        }

        // Convert back to BigInteger
        BigInteger[] result = new BigInteger[n];
        for (int i = 0; i < n; i++)
            result[i] = BigInteger.valueOf((long) Math.round(X[i]));

        return result;
    }
}
