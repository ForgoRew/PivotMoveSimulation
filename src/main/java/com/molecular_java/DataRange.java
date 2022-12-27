package com.molecular_java;

/**
 * A simple class used in the {@code CountErrBlockMethod(SimSpace s)}
 * method to select ranges in data file for each block.
 */
public class DataRange {
    /**Tells in which cycle the block starts. */
    public int from;
    /**Tells in which cycle the block ends. */
    public int to;
    
    /**
     * Straith forward constructor of the data range.
     * From and to cycles are stored in the object.
     * @param From tells from which cycle the block starts
     * @param To tells the number of cycle where the block ends.
     */
    public DataRange(int From, int To){
        from = From;
        to = To;
    }
}
