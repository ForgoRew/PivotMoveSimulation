package com.molecular_java;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
/**
 * This class is used to write to files instead of standard {@link BufferedWriter}. 
 * The {@link BufferedWriter} has not so neat methods. It contains
 * a constructor which gets the initialized {@link BufferedWriter}
 * and loads it as an attribute of this class.
 * 
 * Then it has a {@link #writeLine(String)} and a {@link #writeAll(BufferedReader)}
 * methods which easily write to the loaded buffered writer.
 */
public class MyWriter {
    /** Nested BufferedWriter used to write things. */
    BufferedWriter buffy;
    /**
     * Constructor which takes a {@link BufferedWriter} and
     * uses it as an attribute.
     * @param bufWrite buffered writer used as an attribute
     */
    public MyWriter(BufferedWriter bufWrite){
        buffy = bufWrite;
    }
    /**
     * Writes given String as a line to a {@code buffy}.
     * @param str string to write
     * @throws IOException if an I/O Error occurs
     */
    public void writeLine(String str) throws IOException{
        buffy.write(str);
        buffy.newLine();
    }
    /**
     * Prints all lines of the given file to the {@link #buffy}.
     * @param in given input file
     * @throws IOException if an I/O Error occurs
     */
    public void writeAll(BufferedReader in) throws IOException{
        String line = in.readLine();
        while (line != null){
            this.writeLine(line);
            line = in.readLine();
        }
    }
    /**
     * Closes {@link #buffy}.
     * @throws IOException if an I/O Error occurs
     */
    public void close() throws IOException{
        buffy.close();
    }
}
