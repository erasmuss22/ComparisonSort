import java.util.Random;

/**
 * This program tests some of the functionality of the ComparisonSort class. It
 * does not completely test the ComparisonSort class. You should make sure that
 * you do completely test your ComparisonSort class, either by modifying this
 * file or by writing a different driver.
 */
public class TestSort {

    /**
     * Main method to run the ComparisonSort class.
     * 
     * @param args
     *            a two-value array: first the number of items in the input
     *            array, then the random number seed (integer)to use in
     *            generating values
     */
    public static void main(String[] args) {

        if (args.length != 2) {
            System.err.println("Expected 2 but got " + args.length);
            System.err.println("Arguments expected:");
            System.err.println("  # items in input array");
            System.err.println("  random # seed");
            System.exit(1);
        }
        int arrSize = Integer.parseInt(args[0]);
        int seed = Integer.parseInt(args[1]);

        System.out.format("sorting %,d items generated with random seed %d\n\n", arrSize, seed);

        // Create the input array of unsorted objects.
        SortObject[] arr = new SortObject[arrSize];

        // It is important to give the seed so you can reproduce results.
        Random random = new Random(seed);
        for (int k = 0; k < arrSize; k++)
            arr[k] = new SortObject(random.nextInt());


        // Run all the sorts on the array of random integers.
        ComparisonSort.runAllSorts(arr);

    }
}
