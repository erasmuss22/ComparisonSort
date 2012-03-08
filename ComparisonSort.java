///////////////////////////////////////////////////////////////////////////////
//                   ALL STUDENTS COMPLETE THESE SECTIONS
// Main Class File:  TestSort.java
// File:             ComparisonSort.java
// Semester:         Spring 2011
//
// Author:           Erin Rasmussen  ejrasmussen2@wisc.edu
// CS Login:         rasmusse
// Lecturer's Name:  Beck Hasti
// Lab Section:      Lecture 2
//
//                   
//                   STUDENTS WHO GET HELP FROM ANYONE OTHER THAN THEIR PARTNER
// Credits:          Some code taken from the online readings
//////////////////////////// 80 columns wide //////////////////////////////////
/**
 * This class implements seven different comparison sorts:
 * <ul>
 * <li>selection sort</li>
 * <li>insertion sort</li>
 * <li>merge sort</li>
 * <li>quick sort using median-of-three partitioning</li>
 * <li>quick sort using first-value paritioning</li>
 * <li>heap sort</li>
 * <li>two-way selection sort</li>
 * </ul>
 * It also has a method that runs all the sorts on the same input array and
 * prints out statistics.
 */

public class ComparisonSort {
	//used to store the data moves in a sort
	private static int dataMoves = 0;

	/**
	 * Adds 1 move to dataMoves each time a data move is made
	 *
	 */
	private static void setMoves(){
		dataMoves++;
	}

	/**
	 * Sets the amount of dataMoves to 0 between sorts
	 *
	 */
	private static void resetMoves(){
		dataMoves = 0;
	}

	/**
	 * Returns the amount of dataMoves performed in each sort
	 *
	 * @return an int of the amount of dataMoves
	 */
	private static int getMoves(){
		return dataMoves;
	}
	/**
	 * Sorts the given array using the selection sort algorithm. You may use
	 * either the algorithm discussed in the on-line reading or the algorithm
	 * discussed in lecture (which does fewer data moves than the one from the
	 * on-line reading). Note: after this method finishes the array is in sorted
	 * order.
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 **/
	public static <E extends Comparable<E>> void selectionSort(E[] A) {
		long begin = System.currentTimeMillis();
		int j, k, minIndex;
		E min;
		int N = A.length;
		for (k = 0; k < N; k++) {
			min = A[k];
			setMoves();
			minIndex = k;
			for (j = k+1; j < N; j++) {
				if (A[j].compareTo(min) < 0) {
					min = A[j];
					minIndex = j;
				}
			}
			A[minIndex] = A[k];
			setMoves();
			A[k] = min;
			setMoves();
		}
		long end = System.currentTimeMillis();
		printStatistics("selection", SortObject.getCompares(), getMoves(), 
				end - begin);
		SortObject.resetCompares();
		resetMoves();
	}

	/**
	 * Sorts the given array using the insertion sort algorithm. Note: after
	 * this method finishes the array is in sorted order.
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 **/
	public static <E extends Comparable<E>> void insertionSort(E[] A) {
		long begin = System.currentTimeMillis();
		int k, j;
		E tmp;
		int N = A.length;    
		for (k = 1; k < N; k++) {
			tmp = A[k];
			setMoves();
			j = k - 1;
			while ((j >= 0) && (A[j].compareTo(tmp) > 0)) {
				A[j+1] = A[j]; // move one value over one place to the right
				setMoves();
				j--;
			}
			A[j+1] = tmp;    // insert kth value in correct place relative 
			setMoves();      // to previous values
		}
		long end = System.currentTimeMillis();
		printStatistics("insertion", SortObject.getCompares(), getMoves(), 
				end - begin);
		SortObject.resetCompares();
		resetMoves();

	}

	/**
	 * Sorts the given array using the merge sort algorithm. Note: after this
	 * method finishes the array is in sorted order.
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 **/
	public static <E extends Comparable<E>> void mergeSort(E[] A) {
		long begin = System.currentTimeMillis();
		mergeAux(A, 0, A.length - 1); // call the aux. function to do all
		long end = System.currentTimeMillis();   //the work
		//2 * getMoves() because of arraycopy
		printStatistics("merge", SortObject.getCompares(), 2 * getMoves(),
				end - begin);
		SortObject.resetCompares();
		resetMoves();
	}

	/**
	 * An auxiliary method for mergeSort to allow recursive calls.
	 *
	 * @param A the array to sort
	 * @param low the lowest position in the array
	 * @param high the highest position in the array
	 */
	@SuppressWarnings("unchecked")
	private static <E extends Comparable<E>> void mergeAux(E[] A, int low, 
			int high) {
		// base case
		if (low == high) return;
		// recursive case
		//Step 1: Find the middle of the array (conceptually, divide it in half)
		int mid = (low + high) / 2;
		// Steps 2 and 3: Sort the 2 halves of A
		mergeAux(A, low, mid);
		mergeAux(A, mid+1, high);
		// Step 4: Merge sorted halves into an auxiliary array
		E[] tmp = (E[])(new Comparable[high-low+1]);
		int left = low;    // index into left half
		int right = mid+1; // index into right half
		int pos = 0;       // index into tmp
		while ((left <= mid) && (right <= high)) {
			// choose the smaller of the two values "pointed to" by left, right
			// copy that value into tmp[pos]
			// increment either left or right as appropriate
			// increment pos
			if ((A[left].compareTo(A[right]) < 0)) {
				tmp[pos] = A[left];
				setMoves();
				left++;
			}
			else {
				tmp[pos] = A[right];
				setMoves();
				right++;
			}
			pos++;
		}
		// when one of the two sorted halves has "run out" of values, but
		// there are still some in the other half; copy all the remaining 
		// values to tmp
		// Note: only 1 of the next 2 loops will actually execute
		while (left <= mid) {
			tmp[pos] = A[left]; 
			setMoves();
			left++;
			pos++;
		}
		while (right <= high) {
			tmp[pos] = A[right];
			setMoves();
			right++;
			pos++;
		}      
		System.arraycopy(tmp, 0, A, low, tmp.length);
	}

	/**
	 * Sorts the given array using the quick sort algorithm, using the median of
	 * the first, last, and middle values in each segment of the array as the
	 * pivot value. Note: after this method finishes the array is in sorted
	 * order.
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 **/
	public static <E extends Comparable<E>> void quickSortUsingMedianOfThree
	(E[] A) {
		long begin = System.currentTimeMillis();
		quickAux(A, 0, A.length - 1);
		long end = System.currentTimeMillis();
		printStatistics("quick (median of three)", SortObject.getCompares(), 
				getMoves(), end - begin);
		SortObject.resetCompares();
		resetMoves();
	}

	/**
	 * An auxiliary method for quickSortUsingMedianOfThree to allow recursive 
	 * calls.
	 *
	 * @param A the array to sort
	 * @param low the lowest position in the array
	 * @param high the highest position in the array
	 */
	private static <E extends Comparable<E>> void quickAux(E[] A, int low, 
			int high) {
		if (high - low <= 2 ) {
			int mid = (high + low) / 2;
			if (A[low].compareTo(A[mid]) <= 0 && A[low].compareTo(A[high])
					<= 0){
				if (A[mid].compareTo(A[high]) <= 0)return;
				else {
					swap(A, mid, high);
				}
			}
			else if (A[mid].compareTo(A[low]) <= 0 && A[mid].compareTo(A[high])
					<= 0){
				swap(A, low, mid);
				if (A[mid].compareTo(A[high]) > 0){
					swap(A, mid, high);
				}
			}
			else {
				swap(A, high, low);
				if (A[mid].compareTo(A[high]) > 0){
					swap(A, mid, high);
				}
			}

		}
		else {
			int right = partition(A, low, high);
			quickAux(A, low, right);
			quickAux(A, right + 2, high);
		}
	}

	/**
	 * Partitions the array into parts greater than and less than the pivot
	 * and puts the pivot in its place
	 *
	 * @param A the array to sort
	 * @param low the lowest position in the array
	 * @param high the highest position in the array
	 */
	private static <E extends Comparable<E>> int partition(E[] A, int low, 
			int high) {
		E pivot = pivot(A, low, high);
		int left = low + 1; 
		int right = high - 2;
		while (left <= right) {
			while (A[left].compareTo(pivot) < 0) left++;
			while (A[right].compareTo(pivot) > 0) right--;
			if (left <= right) {
				swap(A, left, right);
				left++;
				right--;
			}
		}
		swap(A, left, high-1);
		return right;
	}

	/**
	 * Swaps the elements between the given positions
	 *
	 * @param A the array to sort
	 * @param val1 the position of the first value
	 * @param val2 the position of the second value
	 */
	private static <E> void swap(E[] A, int val1, int val2){
		E tmp = A[val1];
		setMoves();
		A[val1] = A[val2];
		setMoves();
		A[val2] = tmp;
		setMoves();
	}

	/**
	 * Finds the median of 3 values and puts them and the pivot in the correct
	 * positions
	 *
	 * @param A the array to sort
	 * @param low the lowest position in the array
	 * @param high the highest position in the array
	 * @return the pivot of the array
	 */
	private static <E extends Comparable<E>> E pivot(E[] A, int low, int high){
		int mid = ((high + low) / 2);
		if (A[low].compareTo(A[mid]) <= 0 && A[low].compareTo(A[high]) <= 0){
			if (A[mid].compareTo(A[high]) > 0)swap(A, mid, high);
		}
		else if (A[mid].compareTo(A[low]) < 0 && A[mid].compareTo(A[high]) <=
			0){
			swap(A, low, mid);
			if (A[mid].compareTo(A[high]) > 0)swap(A, mid, high);
		}
		else {
			swap(A, high, low);
			if (A[high].compareTo(A[mid]) < 0) swap(A, mid, high);
		}
		swap(A, mid, high - 1);
		return A[high - 1];
	}


	/**
	 * Sorts the given array using the quick sort algorithm, using the first
	 * value in each segment of the array as the pivot value. Note: after this
	 * method finishes the array is in sorted order.
	 * 
	 * The algorithm for quick sort using the first value as the pivot value is
	 * as follows:
	 * 
	 * <pre>
	 * Partition A[low] ... A[high] using A[low] as the pivot value so that
	 * after partitioning, 
	 *   - A[low]...A[pivot-1] contains values <= the pivot value,
	 *   - A[pivot] contains the pivot value, and
	 *   - A[pivot+1]...A[high] contains values >= the pivot value
	 * 
	 * To partition from low to high:
	 *   (1) Initialize left to low + 1.
	 *   (2) Initialize right to high.
	 *   (3) Initialize the pivot value to the array value in the low slot.
	 * 
	 * Using the same technique done in the median-of-three partitioning,
	 * increment left and decrement right (and swap values as appropriate)
	 * until left and right have crossed.  As you do so, make sure that right
	 * doesn't become less than low.
	 * 
	 * If left becomes greater than high, that means the pivot value is the
	 * largest value, so set right to high.
	 * 
	 * right is now the index where the the pivot value should go, so swap the
	 * pivot value into place and return right.
	 * 
	 * 
	 * To quick sort from low to high:
	 *   (1) Partition between low and high; call the returned index "pivot".
	 *   (2) If the subarray from index low to index pivot-1 contains at least
	 *       one element, quick sort the subarray.
	 *   (3) If the subarray from index pivot+1 to index high contains at least
	 *       one element, quick sort the subarray.
	 * </pre>
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 **/
	public static <E extends Comparable<E>> void quickSortUsingFirst(E[] A) {
		long begin = System.currentTimeMillis();
		quickFirstAux(A, 0, A.length - 1);
		long end = System.currentTimeMillis();
		printStatistics("quick (first)", SortObject.getCompares(), getMoves(),
				end - begin);
		SortObject.resetCompares();
		resetMoves();
	}

	/**
	 * An auxiliary method for quickSortUsingFirst to allow recursive calls.
	 *
	 * @param A the array to sort
	 * @param low the lowest position in the array
	 * @param high the highest position in the array
	 */
	private static <E extends Comparable<E>> void quickFirstAux(E[] A, int low,
			int high) {
		int pivot = partitionFirst(A, low, high);
		if ((pivot - 1 - low) >= 1) quickFirstAux(A, low, pivot - 1);
		if ((high - pivot + 1) > 1)quickFirstAux(A, pivot + 1, high);
	}

	/**
	 * Partitions the array by putting values less than the first position
	 * to the left, and values greater than the pivot to the right and then
	 * places the pivot in between.
	 *
	 * @param A the array to sort
	 * @param low the lowest position in the array
	 * @param high the highest position in the array
	 */
	private static <E extends Comparable<E>> int partitionFirst(E[] A, int low,
			int high) {
		E pivot = A[low];
		int left = low + 1; 
		int right = high;
		while (left <= right) {
			while (A[left].compareTo(pivot) < 0 && left < high) left++;
			while (A[right].compareTo(pivot) > 0 && right > low) right--;
			if (left <= right) {
				swap(A, left, right);
				left++;
				right--;
			}
		}
		if (left > high){
			right = high;
			swap(A, right, low);
			return right;
		}
		else swap(A, low, right);
		return right;
	}

	/**
	 * Sorts the given array using the heap sort algorithm outlined below. Note:
	 * after this method finishes the array is in sorted order.
	 * 
	 * <p>
	 * The heap sort algorithm is:
	 * </p>
	 * 
	 * <pre>
	 * for each i from 1 to the end of the array
	 *     insert A[i] into the heap (contained in A[0]...A[i-1])
	 *     
	 * for each i from the end of the array up to 1
	 *     remove the max element from the heap and put it in A[i]
	 * </pre>
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 **/
	@SuppressWarnings("unchecked")
	public static <E extends Comparable<E>> void heapSort(E[] A) {
		long begin = System.currentTimeMillis();
		E[] tmp = (E[])(new Comparable[A.length + 1]);
		int currentNum = 1;
		for (int i = 0; i < A.length; i++) {
			insertInHeapArray(tmp, A[i], currentNum);
			currentNum++;
		}
		int heapEnd = tmp.length - 1;
		for (int i = A.length - 1; i >= 0; i--){
			A[i] = heapRemoveMax(tmp, heapEnd);
			heapEnd--;
		}
		long end = System.currentTimeMillis();
		printStatistics("heap", SortObject.getCompares(), getMoves(), end -
				begin);
		SortObject.resetCompares();
		resetMoves();
	}

	/**
	 * Inserts an element based on its value into a heap array and makes sure
	 * the heap follows the rules for heaps.
	 *
	 * @param A the heap array
	 * @param data the value to be inserted
	 * @param pos the position of the next open value in the heap
	 */
	public static <E extends Comparable<E>> void insertInHeapArray(E[] A, 
			E data, int pos) {
		A[pos] = data;
		setMoves();
		int parent = pos / 2;
		int current = pos;			//void changes pos
		while (parent > 0 && A[parent].compareTo(A[current]) <= 0) { 
			swap(A, current, parent);
			current = parent;
			parent = current / 2;
		}

	}

	/**
	 * Removes and returns the first value in the heap array and readjusts
	 * the heap using the rules of heaps.
	 *
	 * @param B the the array
	 * @param heapEnd, the value of the last position with an element
	 * @return the maximum value in the heap array
	 */
	public static <E extends Comparable<E>> E  heapRemoveMax(E[] B, int 
			heapEnd) {
		E tmp = B[1];
		setMoves();
		reheapify(B, heapEnd);
		return tmp;
	}

	/**
	 * Adjusts the heap array to follow the rules of heap.
	 *
	 * @param B the heap array
	 * @param heapEnd the last element in the heap array
	 */
	public static <E extends Comparable<E>> void reheapify(E[] B, int 
			heapEnd) {
		B[1] = B[heapEnd];
		B[heapEnd] = null;
		int i = 1;
		boolean done = false;
		while (!done){
			int k = i * 2;
			if (B[k] != null && B[i].compareTo(B[k]) <= 0){
				if ((k + 1) < B.length && B[k + 1] != null && B[k].compareTo
						(B[k + 1]) <= 0){
					swap(B, i, k + 1);
					i = k + 1;
					if ((2 * i) >= B.length)done = true;
				}
				else {
					swap(B, i, k);
					i = k;
					if ((2 * i) >= B.length)done = true;
				}
			}
			else if ((k + 1) < B.length && B[k + 1] != null && B[i].compareTo
					(B[k + 1]) <= 0){
				swap(B, i, k + 1);
				i = k + 1;
				if ((2 * i) >= B.length)done = true;
			}
			else done = true;
		}
	}

	/**
	 * Sorts the given array using the two-way selection sort algorithm outlined
	 * below. Note: after this method finishes the array is in sorted order.
	 * <p>
	 * The two-way selection sort is a bi-directional selection sort that sorts
	 * the array from the two ends towards the center. The two-way selection
	 * sort algorithm is:
	 * </p>
	 * 
	 * <pre>
	 * begin = 0, end = A.length-1
	 * 
	 * // At the beginning of every iteration of this loop, we know that the 
	 * // elements in A are in their final sorted positions from A[0] to A[begin-1]
	 * // and from A[end+1] to the end of A.  That means that A[begin] to A[end] are
	 * // still to be sorted.
	 * do
	 *     use the MinMax algorithm (described below) to find the minimum and maximum 
	 *     values between A[begin] and A[end]
	 *     
	 *     swap the maximum value and A[end]
	 *     swap the minimum value and A[begin]
	 *     
	 *     ++begin, --end
	 * until the middle of the array is reached
	 * </pre>
	 * <p>
	 * The MinMax algorithm allows you to find the minimum and maximum of N
	 * elements in 3N/2 comparisons (instead of 2N comparisons). The way to do
	 * this is to keep the current min and max; then
	 * </p>
	 * <ul>
	 * <li>take two more elements and compare them against each other</li>
	 * <li>compare the current max and the larger of the two elements</li>
	 * <li>compare the current min and the smaller of the two elements</li>
	 * </ul>
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 **/
	public static <E extends Comparable<E>> void twoWaySelectSort(E[] A) {
		long timeBegin = System.currentTimeMillis();
		int begin = 0; 
		int end = A.length - 1;
		int currentMax = A.length - 1;
		int currentMin = 0;
		do {
			int left = begin;
			int right = end;
			while (left <= right) {
				if (A[left].compareTo(A[right]) > 0) {
					if (A[left].compareTo(A[currentMax]) > 0) 
						currentMax = left;
					if (A[right].compareTo(A[currentMin]) < 0) 
						currentMin = right;
				}
				else if(A[left].equals(A[right])){
					if (A[left].compareTo(A[currentMin]) < 0) 
						currentMin = left;
					if (A[right].compareTo(A[currentMax]) > 0) 
						currentMax = right;
				}
				else {
					if (A[right].compareTo(A[currentMax]) > 0) 
						currentMax = right;
					if (A[left].compareTo(A[currentMin]) < 0) 
						currentMin = left;
				}
				left++; 
				right--;
			}
			if (begin != currentMax && end != currentMin){
				swap(A, currentMin, begin);
				swap(A, currentMax, end);
			}
			else if(begin == currentMax && end != currentMin){
				swap(A, currentMax, end);
				swap(A, currentMin, begin);
			}
			else if(begin != currentMax && end == currentMin){
				swap(A, currentMin, begin);
				swap(A, currentMax, end);
			}
			else if(begin == currentMin && end != currentMax){
				swap(A, currentMax, end);
			}
			else if(begin != currentMin && end == currentMax) {
				swap(A, currentMin, begin);
			}
			else if(begin == currentMin && end == currentMax){
				swap(A, currentMin, begin);
				swap(A, currentMax, end);
			}
			else {
				E tmp = A[currentMin];
				setMoves();
				swap(A, currentMax, end);
				A[begin] = tmp;
				setMoves();
			}
			begin++; 
			currentMin = begin;
			end--;
			currentMax = end;
		} while (begin <= end);
		long timeEnd = System.currentTimeMillis();
		printStatistics("2-way selection", SortObject.getCompares(), 
				getMoves(), timeEnd - timeBegin);
		SortObject.resetCompares();
		resetMoves();
	}


	/**
	 * Internal helper for printing rows of the output table.
	 * 
	 * @param sort name of the sorting algorithm
	 * @param compares number of comparisons performed during sort
	 * @param moves number of data moves performed during sort
	 * @param milliseconds time taken to sort, in milliseconds
	 */
	private static void printStatistics(String sort, int compares, int moves,
			long milliseconds) {
		System.out.format("%-23s%,15d%,15d%,15d\n", sort, compares, moves, 
				milliseconds);
	}

	/**
	 * Sorts the given array using the seven different sorting algorithms and
	 * prints out statistics. The sorts performed are:
	 * <ul>
	 * <li>selection sort</li>
	 * <li>insertion sort</li>
	 * <li>merge sort</li>
	 * <li>quick sort using median-of-three partitioning</li>
	 * <li>quick sort using first-value paritioning</li>
	 * <li>heap sort</li>
	 * <li>two-way selection sort</li>
	 * </ul>
	 * <p>
	 * The statistics displayed for each sort are: number of comparisons, number
	 * of data moves, and time (in milliseconds).
	 * </p>
	 * <p>
	 * Note: each sort is given the same array (i.e., in the original order) and
	 * the input array A is not changed by this method.
	 * </p>
	 * 
	 * @param A
	 *            the array to sort
	 **/


	static public void runAllSorts(SortObject[] A) {
		System.out.format("%-23s%15s%15s%15s\n", "algorithm", "data compares", "data moves", "milliseconds");
		System.out.format("%-23s%15s%15s%15s\n", "---------", "-------------", "----------", "------------");
		selectionSort(arrayCopy(A));
		insertionSort(arrayCopy(A));
		mergeSort(arrayCopy(A));
		quickSortUsingMedianOfThree(arrayCopy(A));
		quickSortUsingFirst(arrayCopy(A));
		heapSort(arrayCopy(A));
		twoWaySelectSort(arrayCopy(A));
	}

	/**
	 * Copies all values in the array to be sorted so the array's values
	 * never actually change
	 *
	 * @param A the array to copy
	 * @return a copy of the array
	 */
	@SuppressWarnings("unchecked")
	private static <E extends Comparable<E>> E[] arrayCopy(E[] A) {
		E[] B = (E[]) new Comparable[A.length];
		for (int i = 0; i < A.length; i++) B[i] = A[i];
		return B; 
	}    

}
