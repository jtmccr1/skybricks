package dr.skybricks;

import java.util.Arrays;

public class FixedGridEpochProvider extends EpochProvider.Abstract {
    public static final String FIXED_GRID_EPOCH = "fixedGridEpoch";


    public FixedGridEpochProvider(int numGridPoints, double cutOff) {
        super(FIXED_GRID_EPOCH);
        epochCount = numGridPoints+1;

        gridPoints = new double[numGridPoints+1];

        Arrays.fill(gridPoints, 0);
        for (int pt = 0; pt < numGridPoints; pt++) {
        gridPoints[pt] = (pt+1) * (cutOff / numGridPoints);
        }
        gridPoints[numGridPoints] =Double.POSITIVE_INFINITY;
        epochsKnown=true;
    }

    public FixedGridEpochProvider(Double[] grid){
        super(FIXED_GRID_EPOCH);

        boolean needInf = grid[grid.length-1]!=Double.POSITIVE_INFINITY;

        gridPoints = needInf? new double[grid.length+1]:new double[grid.length] ;
        for(int i = 0; i < grid.length; i++){
            gridPoints[i] = grid[i];
        }
        if(needInf){
            gridPoints[gridPoints.length-1]=Double.POSITIVE_INFINITY;
        }
        epochCount = gridPoints.length;
        epochsKnown=true;
    }

    protected void calculateEpochs(){
        epochsKnown=true;
    };

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState() {

    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState() {

    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void acceptState() {

    }
}
