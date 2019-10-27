package VaR;

import java.io.*;
import java.net.URL;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.HashMap;
/**
 * THE ATTRIBUTES OF THE DATA ARE
 *      ﻿Date,Open,High,Low,Close,Volume
 *      We only use Close
 */
public class getStocks {
    public static BufferedReader getCSVfromURL(String urlstr) throws IOException{
        InputStream is = new URL(urlstr).openStream();
        //https://stackoverflow.com/a/4120954
        BufferedReader csv = new BufferedReader(new InputStreamReader(is, "UTF-8"));
        return csv;
        //https://docs.oracle.com/javase/tutorial/networking/urls/readingURL.html
    }
    public static ArrayList<Double> getStocksfromCSV(BufferedReader in) throws IOException{
        String inputLine;
        ArrayList<Double> alData = new ArrayList<Double>();
        in.readLine();//SKIP HEADER
        while ((inputLine = in.readLine()) != null) {
            String[] strTuple = inputLine.split(",");
            /**Col 0: Date
             * Col 1: Open
             * Col 2: High
             * Col 3: Low
             * Col 4: Close
             * Col 5: Volume
             */
            alData.add(Double.parseDouble(strTuple[4]));
        }
        return alData;
    }
    public static String CalculateYearDiffReturnDateAsString(int intYears){
        DateTimeFormatter formatter = DateTimeFormatter.ofPattern("MMM+dd,+yyyy");
        LocalDateTime yearsAgo = LocalDateTime.now().minusYears(intYears);
        //https://stackoverflow.com/questions/22463062/how-to-parse-format-dates-with-localdatetime-java-8
        String fromstr = yearsAgo.format(formatter);
        //CONVERT COMMA TO UNICODE
        fromstr = fromstr.replace(",","%2C");
        return fromstr;
    }
    public static double[][] main(String[] symbols, int intYears) throws IOException{
        System.out.println("=========================================================================");
        System.out.println("getStocks.java");
        System.out.println("=========================================================================");
        String fromStrAPI = CalculateYearDiffReturnDateAsString(intYears);
        int numSym = symbols.length;
        HashMap<String, ArrayList<Double>> mapStocks = new HashMap<>();
        System.out.println("\tFetching VaR.Historic Stock Data from " + intYears + " years(s) ago:");

        //GET STOCK DATA FOR EACH SYMBOL
        for (int i = 0; i < numSym; i++) {
            System.out.println("\t" + symbols[i]);
            //SET urlStrAPI
            String urlStrAPI = "http://www.google.com/finance/historical?q=" + symbols[i] + "&startdate=" + fromStrAPI + "&output=csv";
            BufferedReader csv = getCSVfromURL(urlStrAPI);
            ArrayList<Double> alData = getStocksfromCSV(csv);
            csv.close();
            int size = alData.size();
            System.out.println("\t\tRetrieved " + size + " rows of Stock data");
            System.out.println("\t\t" + urlStrAPI);
            //CONVERT ARRAY LIST DOUBLE TO ARRAY DOUBLE
            mapStocks.put(symbols[i], alData);
        }
        //LOOP THROUGH EACH SYMBOL AND GET THE MINIMUM SIZE
        //THIS IS IN CASE YOU HAVE MULTIPLE STOCKS AND THE ARRAYLISTS ARE OF DIFFERENT SIZES
        int numTuples = Integer.MAX_VALUE;
        for(int i = 0; i < numSym; i++){
            int checkSize  = new ArrayList<>(mapStocks.get(symbols[i])).size();
            numTuples =  Math.min(checkSize,numTuples);
        }
        double[][] stockPrices = new double[numSym][numTuples];
        //POPULATE stockPrices ARRAY. FEED DATA FROM HashMap mapStocks
        for(int i = 0; i < numSym; i++)
            for(int j =0; j < numTuples; j++) {
                ArrayList<Double> arrayListStockPrices = mapStocks.get(symbols[i]);
                stockPrices[i][j] = arrayListStockPrices.get(j);
            }
        return stockPrices;
    }
}