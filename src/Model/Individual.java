package Model;

/**
 * Created by hp on 2016/2/22.
 */
public class Individual {
    public String geneSerial;
    public double fitness;
    public double[] constants;

    public Individual(int constantSCount) {
        geneSerial="";
        fitness=0;
        constants=new double[constantSCount];
    }
}
