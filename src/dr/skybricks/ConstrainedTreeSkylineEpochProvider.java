package dr.skybricks;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.coalescent.TreeIntervalList;
import dr.evolution.tree.NodeRef;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.XMLParseException;

import java.util.Arrays;

/**
 * A skyline like interval list but the group sizes here are based
 * on subtree roots in a constrained tree. This is for use in
 * large trees with unresolved nodes. The heights of the inserted nodes are not
 * informed by the data outside of the information provided by the clade root.
 *
 *
 */
public class ConstrainedTreeSkylineEpochProvider extends EpochProvider.Abstract {
    private final TreeIntervalList intervalList;
    private final Parameter groupSizeParameter;
    private boolean storedEpochsKnown;
    private double[] storedGridPoints;
    private final ConstrainedTreeModel tree;

    public ConstrainedTreeSkylineEpochProvider(TreeIntervalList intervalList, Parameter groupSizeParameter) throws XMLParseException {
        super(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD); //For fun ;)

      if (!(intervalList.getTree() instanceof ConstrainedTreeModel)) {
          throw new IllegalArgumentException("ConstrainedTreeSkylineEpochProvider requires intervals that wrap a constrained Tree");
      }


        this.tree= (ConstrainedTreeModel) intervalList.getTree();
        if (groupSizeParameter.getDimension() > (intervalList.getSampleCount() - 1)) {
            throw new XMLParseException("There are more groups (" + groupSizeParameter.getDimension()
                    + ") than coalescent nodes in the tree (" + (intervalList.getSampleCount() - 1) + ").");
        }

        this.groupSizeParameter = groupSizeParameter;
        int events  = tree.getSubtreeCount();
        int eventsCovered = 0;
        int groupCount = groupSizeParameter.getDimension();
        for (int i = 0; i < groupCount; i++) {
            eventsCovered += getGroupSize(i);
        }

        if (eventsCovered != events) {

            if (eventsCovered == 0 || eventsCovered == groupCount) {
                double[] uppers = new double[groupCount];
                double[] lowers = new double[groupCount];

                // For these special cases we assume that the XML has not specified initial group sizes
                // or has set all to 1 and we set them here automatically...
                int eventsEach = events / groupCount;
                int eventsExtras = events % groupCount;
                for (int i = 0; i < groupCount; i++) {
                    if (i < eventsExtras) {
                        groupSizeParameter.setParameterValue(i, eventsEach + 1);
                    } else {
                        groupSizeParameter.setParameterValue(i, eventsEach);
                    }
                    uppers[i] = Double.MAX_VALUE;
                    lowers[i] = -Double.MAX_VALUE;
                }

                groupSizeParameter.addBounds(new Parameter.DefaultBounds(uppers, lowers));
            } else {
                // ... otherwise assume the user has made a mistake setting initial group sizes.
                throw new IllegalArgumentException("The sum of the initial group sizes does not match the number of  events in the tree.");
            }
        }


        this.intervalList = intervalList;
        if (intervalList instanceof Model) {
            addModel((Model) intervalList);
        }



        addVariable(groupSizeParameter);


        epochCount = groupSizeParameter.getDimension();
        gridPoints = new double[groupSizeParameter.getDimension()];
        storedGridPoints = new double[groupSizeParameter.getDimension()];

        calculateEpochs();

    }

    protected void calculateEpochs () {
        Arrays.fill(gridPoints, 0.0);
        double timeEnd = 0.0;
        int groupIndex = 0;
        int subIndex = 0;
        for (int i = 0; i < intervalList.getIntervalCount(); i++) {

            timeEnd += intervalList.getInterval(i);

            if (intervalList.getIntervalType(i) == IntervalType.COALESCENT) {
                NodeRef node = intervalList.getCoalescentNode(i);
                TreeModel subtree =  tree.getSubtree(node);
                NodeRef nodeInSubtree = tree.getNodeInSubtree(subtree,node);
                if(subtree.isRoot(nodeInSubtree)){
                    subIndex += 1;
                }
                if (subIndex >= groupSizeParameter.getParameterValue(groupIndex)) {
                    gridPoints[groupIndex] = timeEnd;
                    groupIndex += 1;
                    subIndex = 0;
                }
            }
        }
        epochsKnown = true;
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState () {
        storedEpochsKnown = epochsKnown;
        System.arraycopy(gridPoints, 0, storedGridPoints, 0, gridPoints.length);
    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState () {
        epochsKnown = storedEpochsKnown;
        double[] tmp1 = storedGridPoints;
        storedGridPoints = gridPoints;
        gridPoints = tmp1;
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void acceptState () {

    }
    private   int getGroupSize(int groupIndex) {
        double g = groupSizeParameter.getParameterValue(groupIndex);
        if (g != Math.round(g)) {
            throw new RuntimeException("Group size " + groupIndex + " should be integer but found:" + g);
        }
        return (int)Math.round(g);
    }
}

