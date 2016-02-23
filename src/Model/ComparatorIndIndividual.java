package Model;

import java.util.Comparator;

/**
 * Created by hp on 2016/2/23.
 */
public class ComparatorIndIndividual implements Comparator {
    @Override
    public int compare(Object o1, Object o2) {
        Individual individual1=(Individual)o1;
        Individual individual2=(Individual)o2;
        Double a=(Double)individual1.fitness;
        Double b=(Double)individual2.fitness;
        return b.compareTo(a);
    }
}
