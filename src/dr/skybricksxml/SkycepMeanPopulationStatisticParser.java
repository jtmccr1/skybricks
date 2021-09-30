package dr.skybricksxml;

import dr.skybricks.SkycepLikelihood;
import dr.xml.*;

public class SkycepMeanPopulationStatisticParser extends AbstractXMLObjectParser {
    public final static String SKYCEP_MEAN_POPULATION_STATISTIC = "skycepMeanPopulationStatistic";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SkycepLikelihood skycepLikelihood = (SkycepLikelihood) xo.getChild(SkycepLikelihood.class);
        return skycepLikelihood.getStatistic("populationMean");
    }

    /**
     * @return an array of syntax rules required by this element.
     * Order is not important.
     */
    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SkycepLikelihood.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "returns a statistic that reports the population sizes of skycep model";
    }

    @Override
    public Class getReturnType() {
        return SkycepLikelihood.PopulationMeanStatistic.class;
    }

    /**
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    @Override
    public String getParserName() {
        return SKYCEP_MEAN_POPULATION_STATISTIC;
    }
}


