package dr;

import dr.app.plugin.Plugin;
import dr.skybricks.RootHeightProxyParameter;
import dr.skybricksxml.*;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObjectParser;

import java.util.HashSet;
import java.util.Set;

public class Skybricks implements Plugin {
    @Override
    public Set<XMLObjectParser> getParsers() {
        Set<XMLObjectParser> parsers = new HashSet<XMLObjectParser>();

        AbstractXMLObjectParser epochParser = new EpochParser();
        parsers.add(epochParser);

        AbstractXMLObjectParser epochTimeStatistic = new EpochTimeStatistic();
        parsers.add(epochTimeStatistic);

        AbstractXMLObjectParser exponentialParameterParser = new ExponentialParameterParser();
        parsers.add(exponentialParameterParser);

        AbstractXMLObjectParser gmrfSmoothingPriorParser = new GMRFSmoothingPriorParser();
        parsers.add(gmrfSmoothingPriorParser);

        AbstractXMLObjectParser skybrickLikelihoodParser = new SkybrickLikelihoodParser();
        parsers.add(skybrickLikelihoodParser);

        AbstractXMLObjectParser skycepLikelihoodParser = new SkycepLikelihoodParser();
        parsers.add(skycepLikelihoodParser);

        AbstractXMLObjectParser skycepMeanPopulationStatisticParser = new SkycepMeanPopulationStatisticParser();
        parsers.add(skycepMeanPopulationStatisticParser);

        AbstractXMLObjectParser skycepLogPopulationStatisticParser = new SkycepLogPopulationSizeStatisticParser();
        parsers.add(skycepLogPopulationStatisticParser);

        AbstractXMLObjectParser skyRateParameterParser = new SkyRateParameterParser();
        parsers.add(skyRateParameterParser);



        AbstractXMLObjectParser singleDimStat = new SingleDimensionStatisticParser();
        parsers.add(singleDimStat);

        AbstractXMLObjectParser rootHeightProxyParser =  RootHeightProxyParameter.PARSER;
        parsers.add(rootHeightProxyParser);

        return parsers;
    }
}
