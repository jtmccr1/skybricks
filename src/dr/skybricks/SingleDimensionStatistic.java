package dr.skybricks;

import dr.inference.model.Parameter;
import dr.inference.model.Statistic;

public class SingleDimensionStatistic extends Statistic.Abstract{
    Parameter parameter;
    int dimension;

    public SingleDimensionStatistic(Parameter parameter,int dimension){
        this.parameter = parameter;
        this.dimension = dimension;
    }

    @Override
    public int getDimension() {
        // TODO Auto-generated method stub
        return 1;
    }

    @Override
    public double getStatisticValue(int dim) {
        if(dim!=0){
            throw new IllegalArgumentException("Single Dimension stat has dimension 1");
        }
        return parameter.getParameterValue(this.dimension);
    }
}
