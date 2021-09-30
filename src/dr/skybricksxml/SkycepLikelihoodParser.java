package dr.skybricksxml;

import dr.evolution.coalescent.IntervalList;
import dr.inference.model.Parameter;
import dr.skybricks.EpochProvider;
import dr.skybricks.SkycepLikelihood;
import dr.xml.*;

public class SkycepLikelihoodParser extends AbstractXMLObjectParser {
    public static final String SKYCEP = "skycepLikelihood";
    public static final String INTERVALS = "intervals";


    public static final String FIRST_MU = "firstMu";
    public static final String EPOCHS = "epochs";
    public static final String ALPHA= "alpha";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        XMLObject cxo = xo.getChild(FIRST_MU);
        Parameter firstMu = (Parameter) cxo.getChild(Parameter.class);

//        cxo = xo.getChild(ALPHA);
//        Parameter alpha = (Parameter) cxo.getChild(Parameter.class);
        Parameter alpha = new Parameter.Default(3.0);
        cxo = xo.getChild(EPOCHS);
        EpochProvider epochProvider = (EpochProvider) cxo.getChild(EpochProvider.class);

        IntervalList intervalList = (IntervalList) xo.getChild(INTERVALS).getChild(IntervalList.class);
        return new SkycepLikelihood(epochProvider, intervalList,firstMu,alpha);
    }

    /**
     * @return an array of syntax rules required by this element.
     * Order is not important.
     */
    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(FIRST_MU, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
//                new ElementRule(ALPHA, new XMLSyntaxRule[]{
//                        new ElementRule(Parameter.class)
//                }),
                new ElementRule(EPOCHS, new XMLSyntaxRule[]{
                        new ElementRule(EpochProvider.class)
                }),

                new ElementRule(INTERVALS, new XMLSyntaxRule[]{
                        new ElementRule(IntervalList.class)
                }, "The interval list from which the coalescent likelihood will be calculated")
        };
    }

    @Override
    public String getParserDescription() {
        return "Parses the a flexible piecewise coalescent";
    }

    @Override
    public Class getReturnType() {
        return SkycepLikelihood.class;
    }

    /**
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    @Override
    public String getParserName() {
        return SKYCEP;
    }
}
