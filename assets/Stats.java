package VaR;
import org.apache.commons.math3.linear.*;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.*;

public class Stats {
    //instance variable
    private double[] singleVector;
    private double[] xVector;
    private double[] yVector;
    private double[][] multiVector;
    private int numRow;
    private int numCol;
    //constructor
    public Stats(double[] singleVector){
        this.singleVector = singleVector;
        this.numRow = singleVector.length;
    }
    public Stats(double[] xVector, double[] yVector){
        this.xVector = xVector;
        this.yVector = yVector;
        this.numRow = xVector.length;
    }
    public Stats(double[][] multiVector){
        this.multiVector = multiVector;
        this.numCol = multiVector.length;
        this.numRow = multiVector[0].length;
    }
    //INSTANCE METHODS
    public double getMean(){
        double sum = 0.0;
        for(int i = 0; i < numRow; i++)
            sum += singleVector[i];
        return sum/numRow;
    }

    public double getEWVariance(){
        double sum = 0;
        for (int i = 0; i < numRow; i++)
            sum +=(xVector[i])*(yVector[i]);
        return sum/(numRow -1);
    }
    public double getEWMAVariance(){
        double lambda = 0.94;
        double EWMA = xVector[numRow -1]* yVector[numRow -1];
        for (int i = 1; i < numRow; i++)
            EWMA = lambda * EWMA + (1-lambda) * xVector[numRow -1 - i]* yVector[numRow -1 - i];
        return EWMA;
    }
    public double getGARCH11Variance(){
        double uSquared[] = new double [numRow -1];
        for (int i = 0; i < (numRow -1); i++)
            uSquared[uSquared.length - 1 - i] = xVector[i] * yVector[i];
        double[] parameters = LevenbergMarquardt(uSquared);
        double      omega = parameters[0]
                ,   alpha = parameters[1]
                ,   beta = parameters[2];
        double sigmaSquared = uSquared[0];
        for (int i = 1; i < uSquared.length;i++)
            sigmaSquared = omega + (alpha*uSquared[i]) + (beta*sigmaSquared);
        return sigmaSquared;
    }
    //VOLATILITY MEASURES
    public double getEWVolatility(){return Math.sqrt(getEWVariance());}
    public double getEWMAVolatility(){return Math.sqrt(getEWMAVariance());}
    public double getGARCH11Volatility(){return Math.sqrt(getGARCH11Variance());}

    public double getVariance(int measure){
        if(measure == 1)
            return getEWVariance();
        else if (measure == 2)
            return getEWMAVariance();
        else if (measure == 3)
            return getGARCH11Variance();
        else return 0;
    }
    public double getVolatility(int measure){
        if(measure == 1)
            return getEWVolatility();
        else if (measure == 2)
            return getEWMAVolatility();
        else if (measure == 3)
            return getGARCH11Volatility();
        else return 0;
    }
    public double[][] getCorrelationMatrix(int measure){
        double[][] matrix = new double[numCol][numCol];
            for (int i = 0; i < numCol; i++) {
                for (int j = 0; j < numCol; j++) {
                    double covXY = new Stats(multiVector[i], multiVector[j]).getVariance(measure);
                    double sigmaX = new Stats(multiVector[i], multiVector[i]).getVolatility(measure);
                    double sigmaY = new Stats(multiVector[j], multiVector[j]).getVolatility(measure);
                    matrix[i][j] = covXY / (sigmaX * sigmaY);
                }
                System.out.println("\t\t" + Arrays.toString(matrix[i]));
            }
        return matrix;
    }
    public double[][] getCovarianceMatrix(int measure){
        double[][] covarianceMatrix = new double[numCol][numCol];
            for (int i = 0; i < numCol; i++)
                for (int j = 0; j < numCol; j++)
                    covarianceMatrix[i][j] = new Stats(multiVector[i], multiVector[j]).getVariance(measure);
        return covarianceMatrix;
    }
    public double[][] getCholeskyDecomposition(int measure){
        double[][] covarianceMatrix = getCovarianceMatrix(measure);
        double[][] cholesky = new double[covarianceMatrix.length][covarianceMatrix.length];
        //Start at the top right
        for(int i = 0; i < covarianceMatrix.length; i++)
            for(int j = 0; j <= i; j++) {
                double sum = 0;
                for (int k = 0; k < j; k++)
                    sum += cholesky[i][k] * cholesky[j][k];
                if(i==j)
                    cholesky[i][j] = Math.sqrt(covarianceMatrix[i][j] - sum);
                else
                    cholesky[i][j] = (covarianceMatrix[i][j] - sum) / cholesky[j][j];
            }
        return cholesky;
    }
    public double[][] getPercentageChanges(){
        double[][] priceDiff = new double[numCol][numRow - 1];
        for  (int i = 0; i < numCol; i++)
            for (int j = 0; j < numRow - 1; j++)
                priceDiff[i][j] = ((multiVector[i][j]- multiVector[i][j+1])/ multiVector[i][j+1]);
        return priceDiff;
    }
    public double[][] getAbsoluteChanges(){
        double[][] priceDiff = new double[numCol][numRow - 1];
        for  (int i = 0; i < numCol; i++)
            for (int j = 0; j < numRow - 1; j++)
                priceDiff[i][j] = multiVector[i][j]- multiVector[i][j+1];
        return priceDiff;
    }
    public void printMatrixToCSV(String[] header, String title, String relativePath)throws IOException{
        //https://stackoverflow.com/questions/15364342/export-array-values-to-csv-file-java
        //https://stackoverflow.com/questions/34958829/how-to-save-a-2d-array-into-a-text-file-with-bufferedwriter
        BufferedWriter br = new BufferedWriter(new FileWriter(relativePath + title + ".csv"));
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < header.length ; i++){
            sb.append(header[i]);
            if (i < (header.length - 1))
                sb.append(",");
        }
        sb.append("\n");
        for(int i = 0; i < numRow; i++) {
            for (int j = 0; j < numCol; j++) {
                sb.append(multiVector[j][i] + "");
                if (j < (multiVector.length - 1))
                    sb.append(",");
            }
            sb.append("\n");
        }
        br.write(sb.toString());
        br.close();
    }
    public void printVectorToCSV(String header, String title, String relativePath) throws IOException{
        BufferedWriter br = new BufferedWriter(new FileWriter(relativePath + title + ".csv"));
        StringBuilder sb = new StringBuilder();
        sb.append(header);
        sb.append("\n");
        for(int i = 0; i < numRow; i++) {
            sb.append(singleVector[i]);
            sb.append("\n");
        }
        br.write(sb.toString());
        br.close();
    }
    private int earlyexit = 0;
    private double lambda;
    private double[] LevenbergMarquardt(double[] uSquaredArray){
        lambda = 0.001; //non-dimensional fudge factor
        double[] parameters = new double[3];
        parameters[0] = 0.000001346;//omega
        parameters[1] = 0.08339;//alpha
        parameters[2] = 0.9101;//beta
        double likelihood = likelihood(uSquaredArray,parameters);
        while(true) {
            double[] trialParameters = getTrialParameters(parameters, uSquaredArray);
            double trialLikelihood = likelihood(uSquaredArray, trialParameters);
            if(earlyexit == 1) break;
            if (trialLikelihood > likelihood) {
                parameters = trialParameters;
                lambda *= 0.1; //if successful, use a smaller fudge factor
                likelihood = trialLikelihood;
            }
            else lambda *= 10; //if unsuccessful, use a larger fudge factor
        }
        return parameters;
    }
    private double likelihood(double[] uSquaredArray, double parameters[]){
        double      omega    = parameters[0]
                ,   alpha    = parameters[1] //not to be confused with the alpha matrix in LevenbergMarquardt
                ,   beta     = parameters[2]; //not to be confused with the beta vector in LevenbergMarquardt
        int numTuple = uSquaredArray.length;
        //calculate variance
        double[] variance = new double[numTuple-1];
        variance[0] = uSquaredArray[0];
        for(int i = 1; i < variance.length; i++)
            variance[i] = omega + (alpha * uSquaredArray[i]) + (variance[i-1]* beta);

        double likelihood = 0;
        for(int i = 0; i < variance.length; i++)
            likelihood += -Math.log(variance[i]) - (uSquaredArray[i+1]/variance[i]);
        return likelihood;
    }
    private double[] getTrialParameters(double[] parameters, double[] uSquaredArray){
        double omega    = parameters[0];
        double alpha    = parameters[1]; //not to be confused with the alpha matrix in LevenbergMarquardt
        double beta     = parameters[2]; //not to be confused with the beta vector in LevenbergMarquardt
        int numTuple = uSquaredArray.length;
        //initialize derivatives of likelihood
        double dOmega       = 0.0;
        double dAlpha       = 0.0;
        double dBeta        = 0.0;
        //calculate variance
        double[] variance = new double[numTuple-1];
        variance[0] = uSquaredArray[0];
        for(int i = 1; i < variance.length; i++)
            variance[i] = omega + (alpha * uSquaredArray[i]) + (variance[i-1] * beta);
        //calculate derivatives
        for(int i = 1; i < variance.length; i++) {
            dOmega  += ((-1/variance[i])               + (uSquaredArray[i]/Math.pow(variance[i],2)));
            dAlpha  += (-uSquaredArray[i]/variance[i]) + (Math.pow(uSquaredArray[i],2)/Math.pow(variance[i],2));
            dBeta   += (-variance[i-1]/variance[i])    + ((uSquaredArray[i]*variance[i-1])/Math.pow(variance[i],2));
        }
        //populate vector Beta
        double[] vectorBeta = {-0.5*dOmega, 0.5*dAlpha, -0.5*dBeta};
        //System.out.println(Arrays.toString(vectorBeta));
        double[] trialParameters = new double[parameters.length];
        while(true) {
            //populate curvature matrix
            double[][] curvatureMatrix = {
                        {0.5*dOmega*dOmega * (1 + lambda),   0.5*dOmega*dAlpha,              0.5*dOmega*dBeta}
                    ,   {0.5*dAlpha*dOmega,             0.5*dAlpha*dAlpha * (1 + lambda),    0.5*dAlpha*dBeta}
                    ,   {0.5*dBeta*dOmega,              0.5*dBeta*dAlpha,                    0.5*dBeta*dBeta * (1 + lambda)}
            };
            //solve simultaneous equations to calculate deltaParameters
            RealMatrix coefficients = new Array2DRowRealMatrix(curvatureMatrix);
            DecompositionSolver solver = new LUDecomposition(coefficients).getSolver();
            RealVector constants = new ArrayRealVector(vectorBeta);
            RealVector solution = solver.solve(constants);
            //evaluate χ2(a + δa)
            //Stopping rule. If parameter changes are << 1 then stop.
            if ((solution.getEntry(0) + solution.getEntry(1) + solution.getEntry(2) ) < 0.00001){
                earlyexit = 1;
                break;
            }
            for (int i = 0; i < parameters.length; i++)
                trialParameters[i] = parameters[i] + solution.getEntry(i);
            if (trialParameters[0] < 0.0 //omega less than 0
                    || trialParameters[1] < 0.0 //alpha less than zero
                    || trialParameters[1] > 1.0 //alpha greater than 1
                    || trialParameters[2] < 0.0 //beta less than zero
                    || trialParameters[2] > (1 - trialParameters[1])
                    || trialParameters[1] + trialParameters[2] > 1
                    ) {
                lambda *= 10; //use a larger fudge factor
                continue;
            }
            else
                break;
        }
        return trialParameters;
    }




}

