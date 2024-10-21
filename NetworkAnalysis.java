import java.util.Arrays;
import java.util.Collections;

public class NetworkAnalysis {
	private static int n = 0;
	private static int m = 0;
	private static int h;
	private final static String LIST = "-> ";

	public static void main(String[] args) {
		Integer[] sequences = {2, 5, 3, 2, 7, 2, 4, 3}; // INSERT DEGREES HERE
		if (args.length > 0) {
			sequences = new Integer[args.length];
			for (int i = 0; i < args.length; i++) {
				sequences[i] = Integer.parseInt(args[i]);
			}
		}
		
		int[] torus = {0,1,0}; // INSERT SMALL-WORLD (n,1,r) , if you dont care please insert {0,0,0}
		
		
		Arrays.sort(sequences, Collections.reverseOrder());
		System.out.println("Your sequence: " + Arrays.toString(sequences) + "");
		n = sequences.length;
		m = 0;
		for (int i = 0; i < n; i++) {
			m += sequences[i];
		}
		m /= 2;
		System.out.println("Number of nodes: " + n);
		System.out.println("Number of edges: " + m + "\n");
		
		
		
		// check if sequence is graphic
		if (!isGraphicSequence(Arrays.copyOf(sequences, sequences.length))) {
			System.out.println(LIST + "Sequence is NOT graphic ❌");
			return; // Sequence stellt keinen Graph dar, macht keinen Sinn weiter zu machen
		} else
			System.out.println(LIST + "Sequence is graphic ✅");
		System.out.println();

		
		// durfee number / h-index
		h = durfee(Arrays.copyOf(sequences, sequences.length));
		System.out.println(LIST + "Durfee number/h-index:\th = " + h);
		System.out.println();

		
		System.out.println(LIST + "Check if erdosGallai holds for k = 1, 2, ..., n = " + n);
		erdosGallai(sequences);
		System.out.println("(check yourself if graph is threshold graph: first h inequalities have to be equal, the next have to be <=)");
		System.out.println();


		// check if graph is splitgraph
		if(splitGraph(Arrays.copyOf(sequences, sequences.length), durfee(Arrays.copyOf(sequences, sequences.length))))
			System.out.println(LIST + "Graph is a split Graph ✅");
		else
			System.out.println(LIST + "Graph is NOT a split Graph ❌");
		System.out.println();
		
		System.out.println(LIST + "Splittance number: " + (splittanceNumber(Arrays.copyOf(sequences, sequences.length), durfee(Arrays.copyOf(sequences, sequences.length)))));
		System.out.println();
		
		
	
		// Random Graph Model density
		double degMedian = (2 * m / (double)n);
		System.out.println(LIST + "(Random Graph Model) Density of Graph: " + degMedian / (n - 1));




		System.out.println("\n=== MORE THINGS ===\n");
		
		System.out.print("<deg> = " + degMedian + "\n");
		
		if(torus[2] <= (torus[0]-1)/(4)) {
			System.out.println("Torus Lemma can be used :)");
			
			System.out.println("\t=> check condition <deg> = 2r: " + (2 * m / n) + " =? " + 2*torus[2]);
			System.out.println("\t=> check condition <cc> = (3r-3)/(4r-2): " + "'by hand'" + " =? " + (3*torus[2]-3)/(4*torus[2]-2));
			
			
		}else{
			System.out.println("Torus Lemma CANNOT be used :(");
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


	private static void erdosGallai(Integer[] sequence) {
		// sequence IS ASSUMED TO BE GRAPHIC AND SORTED IN DESCENDING ORDER
		// for correct threshold graph calculation, this function has to be called after the global h has been set correctly

		int[] d = new int[sequence.length];
		for (int i = 0; i < sequence.length; i++) {
			d[i] = sequence[i];
		}

		boolean thresholdGraph = true;
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

			System.out.print("\tk = " + k + ":\t");
			if (leftSum > rightSum) {
				System.out.println("The sequence is not graphic at k = " + k + " ❌");
				thresholdGraph = false;
			} else {
				if(leftSum == rightSum) {
					System.out.print(leftSum + " == " + rightSum + " ✅ (equality)");
					if (k > h) thresholdGraph = false;
				} else {
					System.out.print(leftSum + " <= " + rightSum + " ✅");
					if (k <= h) thresholdGraph = false;
				}
				System.out.println();
			}
		}
		if (thresholdGraph) System.out.println("\n" + LIST + "Graph is a threshold graph ✅");
		else System.out.println("\n" + LIST + "Graph is NOT a threshold graph ❌");
	}

	private static boolean isGraphicSequence(Integer[] sequence) {
		// sequence IS ASSUMED TO BE SORTED IN DESCENDING ORDER
		// Havel Hakimi: https://www.geogebra.org/m/ekajwspy
		
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
