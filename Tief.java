public class Tief {
	public static void main(String[] args) {
		if (args.length != 2) {
			System.out.println("<< This program expects 2 arguments of type int >>");
			System.out.println("args[0] = n, args[1] = k\t==> result = n choose k");
			return;
		}
		int n = Integer.parseInt(args[0]);
		int k = Integer.parseInt(args[1]);

		int result = factorial(n) / (factorial(k) * factorial(n - k));
		System.out.println("n = " + n + ", k = " + k);
		System.out.println("(n choose k) = " + result);
	}

	private static int factorial(int number) {
		int result = 1;

		for (int factor = 2; factor <= number; factor++) {
			result *= factor;
		}

		return result;
	}
}
