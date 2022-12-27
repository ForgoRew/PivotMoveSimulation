package com.molecular_java;

/**
 * A class for time meassurement.
 * It has methods for start and stop of the time
 * period.
 */
public class MyWatches {
    /**Starting time in millis. Set in the {@link #start()} method. */
    double startTime;
    /**Ending time in millis. Set in the {@link #end()} method. */
    double endTime;
    /**Elapsed time in millis (count in the {@link #end()}) method. */
    public double timeElapsed;

    /**
     * On start, the method notes the current time for further counting.
     */
    public void start(){
        startTime = System.currentTimeMillis();
    }
    /**
     * The method notes the current time as the end time and sets the
     * elapsed time as a difference of end and start time.
     */
    public void end(){
        endTime = System.currentTimeMillis();
        timeElapsed = endTime - startTime;
    }
}
