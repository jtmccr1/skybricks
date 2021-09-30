package dr.skybricksxml;

import dr.skybricks.SkycepLikelihood;
import dr.xml.*;

public class SkycepLogPopulationSizeStatisticParser extends AbstractXMLObjectParser {
    public final static String SKYCEP_LOG_POPULATION_SIZE_STATISTIC = "skycepLogPopulationSizeStatistic";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SkycepLikelihood skycepLikelihood = (SkycepLikelihood) xo.getChild(SkycepLikelihood.class);
        return skycepLikelihood.getStatistic("logPopsize");
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
        return "returns a statistic that reports the log population sizes of skycep model";
    }

    @Override
    public Class getReturnType() {
        return SkycepLikelihood.LogPopulationSizeStatistic.class;
    }

    /**
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    @Override
    public String getParserName() {
        return SKYCEP_LOG_POPULATION_SIZE_STATISTIC;
    }
}


