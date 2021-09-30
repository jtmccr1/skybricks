package dr.skybricks;

import dr.evolution.coalescent.IntervalList;
import dr.evomodel.coalescent.AbstractCoalescentLikelihood;
import dr.evomodel.coalescent.CoalescentIntervalProvider;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.inference.model.Variable;
import dr.math.MathUtils;
import dr.skybricksxml.SkycepLikelihoodParser;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.XMLParseException;

import java.util.Arrays;
import java.util.List;

public class SkycepLikelihood extends AbstractCoalescentLikelihood
        implements CoalescentIntervalProvider, Citable {
    private final Parameter firstMu;
    private final EpochProvider epochProvider;
    private final IntervalList intervalsList;

    private final Parameter alpha;

    private final double[] Q;
    private final double[] R;
    private final double[] beta;
    private final double[] mu;
    private final double[] numCoalEvents;

    protected int numEpochs;

    public SkycepLikelihood(EpochProvider epochProvider, IntervalList intervals,Parameter firstMu,Parameter alpha) throws XMLParseException {
        super(SkycepLikelihoodParser.SKYCEP,intervals);

        this.intervalsList = intervals;

        this.epochProvider = epochProvider;
        if(epochProvider instanceof Model){
            addModel((Model) epochProvider);
        }
        numEpochs = this.epochProvider.getEpochCount();
        this.numCoalEvents = new double[numEpochs];
        this.firstMu = firstMu;
        this.alpha = alpha;
        this.Q = new double[numEpochs];
        this.R = new double[numEpochs];
        this.beta = new double[numEpochs];
        this.mu = new double[numEpochs];

        addVariable(firstMu);
        addVariable(alpha);
        calculateLogLikelihood();

        addStatistic(new PopulationMeanStatistic());
        addStatistic(new LogPopulationSizeStatistic());

    }



    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     */
    @Override
    protected double calculateLogLikelihood() {
//        double coalescentLikelihood = piecewiseConstant ? calculateConstantCoalescentLogLikelihood() : calculateExpCoalescentLogLikelihood();
        // Matrix operations taken from block update sampler to calculate data likelihood and field prior
        setupSufficientStatistics();


        double currentLike = 0;
        double a = alpha.getParameterValue(0);
        for (int i = 0; i < numEpochs; i++) {
            currentLike += -(a + numCoalEvents[i]) * Math.log(beta[i] + R[i]) + a * Math.log(beta[i]) - numCoalEvents[i] + Q[i];
        }

        assert capturedAllEvents();

        return currentLike;

    }
    private boolean capturedAllEvents(){
        int captured =0;
        for (int i = 0; i < numCoalEvents.length; i++) {
            captured+=numCoalEvents[i];
        }
        return captured == getNumberOfCoalescentEvents();
    }

    private void setupBetaAndMu(){
        double a=alpha.getParameterValue(0);
        mu[0] = firstMu.getParameterValue(0);
        beta[0] = firstMu.getParameterValue(0)*(a-1);

        for (int i = 1; i < numEpochs; i++) {
            beta[i] = mu[i-1]*(a-1);
            mu[i] = (beta[i]+R[i])/(a+numCoalEvents[i]-1);
        }
    }

    protected void setupSufficientStatistics() {

        Arrays.fill(numCoalEvents, 0.0);
        Arrays.fill(Q, 0.0);
        Arrays.fill(R, 0.0);
//        Arrays.fill(ploidySums, 0);
        //index of smallest grid point greater than at least one sampling/coalescent time in current tree
        int minEpoch;
        //index of greatest grid point less than at least one sampling/coalescent time in current tree
        int maxEpoch;
        int numLineages;
        int currentEpoch;
        int currentTimeIndex;

        double[] currentAndNextTime = new double[2];
        double a = alpha.getParameterValue(0);

        //time of last coalescent event in tree
        double lastCoalescentTime;

        currentTimeIndex = moveToNextTimeIndex(0, currentAndNextTime);

        //  numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
        numLineages = intervalsList.getLineageCount(currentTimeIndex);
        minEpoch = epochProvider.getEpoch(currentAndNextTime[0]);

        currentEpoch = minEpoch;

        lastCoalescentTime = currentAndNextTime[0] + intervalsList.getTotalDuration();

//            theLastTime = lastCoalescentTime;

        maxEpoch = epochProvider.getEpoch(lastCoalescentTime);


        double currentEpochEndTime = epochProvider.getEpochEndTime(currentEpoch);
        double currentEpochStartTime = epochProvider.getEpochStartTime(currentEpoch);

        if (currentEpoch < maxEpoch) {

            while (currentAndNextTime[1] <= currentEpochEndTime) { // epochs are inclusive

                //check to see if interval ends with coalescent event
                //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                    numCoalEvents[currentEpoch] += 1;
                    Q[currentEpoch] += Math.log(a+1);
                }
                R[currentEpoch] +=  (currentAndNextTime[1]- currentAndNextTime[0]) * numLineages * (numLineages - 1) * 0.5;
                currentTimeIndex++;
                currentTimeIndex = moveToNextTimeIndex(currentTimeIndex, currentAndNextTime);
                numLineages = intervalsList.getLineageCount(currentTimeIndex);
            }
            R[currentEpoch] += (currentEpochEndTime-currentAndNextTime[0]) * numLineages * (numLineages - 1) * 0.5;
            currentEpoch++;


            while (currentEpoch < maxEpoch ) {

                //from likelihood of interval between first sampling time and minEpoch

                currentEpochEndTime = epochProvider.getEpochEndTime(currentEpoch);
                currentEpochStartTime = epochProvider.getEpochStartTime(currentEpoch);

                if (currentAndNextTime[1] > currentEpochEndTime) {
                    R[currentEpoch] += (currentEpochEndTime - currentEpochStartTime ) * numLineages * (numLineages - 1) * 0.5;
                    currentEpoch++;

                } else {
                    //check to see if interval ends with coalescent event
                    //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                    if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                        numCoalEvents[currentEpoch] += 1;
                        Q[currentEpoch] += Math.log(a+1);
                    }
                    R[currentEpoch] +=  (currentAndNextTime[1]-currentEpochStartTime) * numLineages * (numLineages - 1) * 0.5;
                    currentTimeIndex++;
                    currentTimeIndex = moveToNextTimeIndex(currentTimeIndex, currentAndNextTime);
                    // numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
                    numLineages = intervalsList.getLineageCount(currentTimeIndex);

                    while (currentAndNextTime[1] <= currentEpochEndTime && currentTimeIndex < intervalsList.getIntervalCount()) {
                        //check to see if interval ends with coalescent event
                        //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                        if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                            numCoalEvents[currentEpoch] += 1;
                            Q[currentEpoch] += Math.log(a+1);

                        }
                        R[currentEpoch] +=  (currentAndNextTime[1]- currentAndNextTime[0]) * numLineages * (numLineages - 1) * 0.5;
                        currentTimeIndex++;
                        currentTimeIndex = moveToNextTimeIndex(currentTimeIndex, currentAndNextTime);
                        numLineages = intervalsList.getLineageCount(currentTimeIndex);
                    }
                    R[currentEpoch] +=  (currentEpochEndTime - currentAndNextTime[0]) * numLineages * (numLineages - 1) * 0.5;
                    currentEpoch++;

                }

            }
            currentEpochStartTime = epochProvider.getEpochStartTime(currentEpoch);
            if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                numCoalEvents[currentEpoch] += 1;
                Q[currentEpoch] += Math.log(a+1);
            }
            R[currentEpoch] +=  (currentAndNextTime[1] - currentEpochStartTime) * numLineages * (numLineages - 1) * 0.5;

            currentTimeIndex++;
            while ((currentTimeIndex) < intervalsList.getIntervalCount()) {

                currentTimeIndex = moveToNextTimeIndex( currentTimeIndex, currentAndNextTime);

                //numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
                numLineages = intervalsList.getLineageCount(currentTimeIndex);
                //check to see if interval is coalescent interval or sampling interval

                //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                    numCoalEvents[currentEpoch] += 1;
                    Q[currentEpoch] += Math.log(a+1);
                }
                R[currentEpoch] +=  (currentAndNextTime[1]- currentAndNextTime[0]) * numLineages * (numLineages - 1) * 0.5;
                currentTimeIndex++;

            }

        }else{
            // only one epoch
            while (currentTimeIndex < intervalsList.getIntervalCount()) {
                //check to see if interval ends with coalescent event
                //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                    numCoalEvents[currentEpoch] += 1;
                    Q[currentEpoch] += Math.log(a+1);
                }
                R[currentEpoch] +=  (currentAndNextTime[1]- currentAndNextTime[0]) * numLineages * (numLineages - 1) * 0.5;
                currentTimeIndex++;
                if ((currentTimeIndex) < intervalsList.getIntervalCount()) {
                    currentTimeIndex = moveToNextTimeIndex( currentTimeIndex, currentAndNextTime);

                    // numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
                    numLineages = intervalsList.getLineageCount(currentTimeIndex);
                }
            }

        }
        setupBetaAndMu();
    }

    //Overwrite could move to abstract coalescent likelihood in future
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    }


    private int moveToNextTimeIndex( int lastTimeIndex, double[] times) {
        int currentTimeIndex = lastTimeIndex;
        double currentTime = intervalsList.getIntervalTime(currentTimeIndex);
        double nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
        while (nextTime <= currentTime && currentTimeIndex + 2 < intervalsList.getIntervalCount()) {
            currentTimeIndex++;
            currentTime = intervalsList.getIntervalTime(currentTimeIndex);
            nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
        }
        times[0] = currentTime;
        times[1] = nextTime;
        return currentTimeIndex;
    }


    @Override
    public int getNumberOfCoalescentEvents() {
        return intervalsList.getIntervalCount() - intervalsList.getSampleCount() + 1;
    }

    @Override
    public double getCoalescentEventsStatisticValue(int i) {
        throw new RuntimeException("What should this be?");
    }


    //Citable implementation
    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "BICEP based coalescent generalized over other epoch models";
    }

    @Override
    public List<Citation> getCitations() {
        return Arrays.asList(new Citation(
                new Author[]{
                        new Author("Remco R", "Bouckaert"),
                },
                "An efficient coalescent epoch model for Bayesian phylogenetic inference",
                2021,
                "BioRxiv",
                67, 719, 728,
                " https://doi.org/10.1101/2021.06.28.450225"
        ));
    }

    /**
     * @return the units for this object.
     */
    @Override
    public Type getUnits() {
        return null;
    }

    /**
     * Sets the units for this object.
     *
     * @param units to use
     */
    @Override
    public void setUnits(Type units) {

    }


    public class PopulationMeanStatistic extends Statistic.Abstract {

        public PopulationMeanStatistic() {
            super("populationMean");
        }

        public int getDimension() { return numEpochs; }

        public double getStatisticValue(int i) {
            return mu[i];
        }
    }

    public class LogPopulationSizeStatistic extends Statistic.Abstract{
        public LogPopulationSizeStatistic() {
            super("logPopsize");
        }

        public int getDimension() { return numEpochs; }

        public double getStatisticValue(int i){
            double sample = MathUtils.nextGamma(alpha.getParameterValue(0)+numCoalEvents[i],beta[i]+R[i]);
            //TODO sample that value
            return Math.log(1/sample);
        }

    }
}


