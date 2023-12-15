package dr.skybricksxml;

import dr.inference.model.Parameter;
import dr.skybricks.SingleDimensionStatistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.AttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public  class SingleDimensionStatisticParser extends AbstractXMLObjectParser {
    public final static String SINGLEDIMSTATISTIC_STRING = "singleDimensionStatistic";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
       Parameter p = (Parameter) xo.getChild(Parameter.class);
         int dimension = xo.getIntegerAttribute("dimension");
         return new SingleDimensionStatistic(p,dimension);
    }

    /**
     * @return an array of syntax rules required by this element.
     * Order is not important.
     */
    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(Parameter.class),
                 AttributeRule.newIntegerRule("dimension")
        };
    }

    @Override
    public String getParserDescription() {
        return "returns a statistic that reports the transition times of an epoch model";
    }

    @Override
    public Class getReturnType() {
        return SingleDimensionStatistic.class;
    }

    /**
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    @Override
    public String getParserName() {
        return SINGLEDIMSTATISTIC_STRING;
    }
}

