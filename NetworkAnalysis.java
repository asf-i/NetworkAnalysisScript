import java.util.Arrays;
import java.util.Collections;

public class NetworkAnalysis {

	public static void main(String[] args) {
		final String LIST = "-> ";
		Integer[] sequences = {2, 5, 3, 2, 7, 2, 4, 3};
		//Integer[] sequences = {2,3,2,7,4,5,3,2,4,2}; // INSERT THE DEGREES 
		
		int[] torus = {0,1,0}; // INSERT SMALL-WORLD (n,1,r) , if you dont care please insert {0,0,0}
		
		
		Arrays.sort(sequences, Collections.reverseOrder());
		System.out.println("Your sequence: " + Arrays.toString(sequences) + "");
		int n = sequences.length;
		int m = 0;
		for (int i = 0; i < n; i++) {
			m += sequences[i];
		}
		m /= 2;
		System.out.println("Number of nodes: " + n);
		System.out.println("Number of edges: " + m + "\n");
		
		
		
		// check if sequence is graphic
		if (isGraphicSequence(Arrays.copyOf(sequences, sequences.length)))
			System.out.println(LIST + "Sequence is graphic ✅");
		else
			System.out.println(LIST + "Sequence is NOT graphic ❌");
		System.out.println();

		
		// durfee number / h-index
		System.out.println(LIST + "Durfee number = h-index = " + durfee(Arrays.copyOf(sequences, sequences.length)));
		System.out.println();

		
		System.out.println("check if erdosGallai holds for k");
		erdosGallai(sequences);
		System.out.println();

		System.out.println("check if Graph is a splitGraph...");
		if(splitGraph(Arrays.copyOf(sequences, sequences.length), durfee(Arrays.copyOf(sequences, sequences.length))))
			System.out.println("splittable ✅");
		else
			System.out.println("NOT splittable ❌");
		System.out.println();
		
		
		System.out.println("The splittance number: " + (splittanceNumber(Arrays.copyOf(sequences, sequences.length), durfee(Arrays.copyOf(sequences, sequences.length)))));
		System.out.println();
		
		System.out.println("check yourself if graph is threshold graph: first h inequalities have to be equal, the next have to be <=");
	
		System.out.println();
		System.out.println("=== MORE THINGS ===");
		System.out.println();
		
		System.out.print("<deg> = " + (2 * m / n) + "\n");
		
		if(torus[2] <= (torus[0]-1)/(4)) {
			System.out.println("Torus Lemma can be used :)");
			System.out.println("=>");
			
			System.out.println("check condition <deg> = 2r: " + (2 * m / n) + " =? " + 2*torus[2]);
			System.out.println("check condition <cc> = (3r-3)/(4r-2): " + "'by hand'" + " =? " + (3*torus[2]-3)/(4*torus[2]-2));
			
			
		}else{
			System.out.println("Torous Lemma CANNOT be used :(");
		}
	
	}
	

	
	public static double splittanceNumber(Integer[] sequence, int h) {
	    // Convert Integer[] to int[] if needed
	    int[] d = new int[sequence.length];
	    for (int i = 0; i < sequence.length; i++) {
	        d[i] = sequence[i];
	    }

	    // Left-hand side: sum of the first h elements
	    int leftSum = 0;
	    for (int i = 0; i < h; i++) {
	        leftSum += d[i];
	    }

	    // Right-hand side: h * (h - 1) + sum of elements from h+1 to n
	    int rightSum = h * (h - 1);
	    for (int j = h; j < d.length; j++) {
	        rightSum += d[j];
	    }

	    // Calculate the splittance using the correct formula
	    double splittance = 0.5 * (rightSum - leftSum);

	    return splittance;
	}
	
	public static boolean splitGraph(Integer[] sequence, int h) {
		
		int[] d = new int[sequence.length];
		for (int i = 0; i < sequence.length; i++) {
			d[i] = sequence[i];
		}
		
		// Left-hand side: sum of the first h elements
        int leftSum = 0;
        for (int i = 0; i < h; i++) {
            leftSum += d[i];
        }

        // Right-hand side: h * (h - 1) + sum of elements from h+1 to n
        int rightSum = h * (h - 1);
        for (int j = h; j < d.length; j++) {
            rightSum += d[j];
        }

        // Check if both sides are equal
        return leftSum == rightSum;
		
	}

	public static int durfee(Integer[] sequence) {

		int[] d = new int[sequence.length];
		for (int i = 0; i < sequence.length; i++) {
			d[i] = sequence[i];
		}

		int h = -1; // Initialize h to -1 (default invalid value)

		// Iterate through the array to find the max index i where d[i] >= i - 1
		for (int i = 1; i <= d.length; i++) {
			if (d[i - 1] >= i - 1) {
				h = Math.max(i, h); // Update h to the current index
			}
		}

		return h;
	}


	public static void erdosGallai(Integer[] sequence) {

		int[] d = new int[sequence.length];
		for (int i = 0; i < sequence.length; i++) {
			d[i] = sequence[i];
		}

		int n = d.length;

		// Sort the sequence in non-increasing order
		Arrays.sort(d);
		for (int i = 0; i < n / 2; i++) {
			int temp = d[i];
			d[i] = d[n - i - 1];
			d[n - i - 1] = temp;
		}

		// Check if the sum is even
		int sum = Arrays.stream(d).sum();
		if (sum % 2 != 0) {
			System.out.println("-> The sequence does not have an even sum => CANNOT BE GRAPHIC");
			System.out.println();
		}

		// Check the prefix sum condition for all k from 1 to n-1
		for (int k = 1; k <= n; k++) {
			int leftSum = 0;
			for (int i = 0; i < k; i++) {
				leftSum += d[i];
			}

			int rightSum = k * (k - 1);
			for (int j = k; j < n; j++) {
				rightSum += Math.min(k, d[j]);
			}
			
			System.out.print(leftSum + " ? " + rightSum + "  <=>  ");

			if (leftSum > rightSum) {
				System.out.println("The sequence is not graphic at k = " + k + " ❌");
			} else {
				System.out.print("holds for k = " + k + " ✅");
				if(leftSum == rightSum) {
					System.out.print(" -> EVEN WITH EQUALITY!");
				}
				System.out.println();
			}
		}
	}

	public static boolean isGraphicSequence(Integer[] sequence) {
		// Havel Hakimi: https://www.geogebra.org/m/ekajwspy
		
		// Sort the sequence in non-increasing order
		Arrays.sort(sequence, Collections.reverseOrder());

		// Continue the algorithm until all degrees are processed
		while (true) {
			// Remove the largest degree (first element after sorting)
			if (sequence[0] == 0) {
				return true; // All remaining degrees are zero, so it's graphic
			}

			int maxDegree = sequence[0];
			if (maxDegree < 0 || maxDegree >= sequence.length) {
				return false; // Invalid sequence if maxDegree is negative or larger than the number of
								// remaining nodes
			}

			// Subtract the degree from the remaining elements
			for (int i = 1; i <= maxDegree; i++) {
				sequence[i]--;
				if (sequence[i] < 0) {
					return false; // Sequence is not graphic if any degree becomes negative
				}
			}

			// Set the first degree to zero and sort the sequence again
			sequence[0] = 0;
			Arrays.sort(sequence, Collections.reverseOrder());
		}

	}

}
