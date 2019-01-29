package omegaPlots;

import java.awt.Color;
import java.util.ArrayList;

public class Calculator {
    static double c = 3e8;
   
    double omega(double k,double omegaP,double epsilon) {
        double c2 = c*c;
        double k2 = k*k;
        double wP2 = omegaP*omegaP;
        double term1 = 0.5*wP2;
        double term2 = k2*c2*((epsilon+1)/(2*epsilon));
        double val = (epsilon*wP2+k2*c2*(epsilon+1));
        double term3 = -(1/(2*epsilon))*Math.sqrt(val*val - 4*epsilon*wP2*k2*c2); 
        return Math.sqrt(term1 + term2 + term3);
    }

    class Equation implements Function {
        double omegaP, k, elow, ehi, d;

        public void setValues(double omegaP,double k, double elow, double ehi,double d) {
            this.omegaP = omegaP;
            this.k = k;
            this.elow = elow;
            this.ehi = ehi;
            this.d = 2*Math.PI*d*c/omegaP;
        }
        
        @Override
        public Complex eval(Complex omegaB) {
            Complex kB = Complex.mul(omegaB,1/c);
            Complex kB2 = Complex.mul(kB,kB);
            Complex z2 = new Complex(k*k);
            Complex wP2 = new Complex(omegaP*omegaP);
            Complex wB2 = Complex.mul(omegaB,omegaB);
            Complex ep = Complex.add(Complex.mul(Complex.div(wP2,wB2),-1),1);
            Complex ulow = Complex.mul(Complex.sqrt(Complex.add(z2,Complex.mul(kB2,-elow))),1/elow);
            Complex uhi = Complex.mul(Complex.sqrt(Complex.add(z2,Complex.mul(kB2,-ehi))),1/ehi);
            Complex up = Complex.div(Complex.sqrt(Complex.sub(z2,Complex.mul(kB2,ep))),ep);
            
            Complex term1 = Complex.tanh(Complex.mul(uhi,ehi*d));
            Complex num = Complex.add(Complex.div(up,ulow),1);
            Complex den = Complex.add(Complex.div(uhi,ulow),Complex.div(up,uhi));
            Complex term2 = Complex.div(num,den);
            return Complex.add(term1, term2);
        }
        
    }

    public void plot(OmegaPlots frame) {
        
        Plot plt=new Plot();
        plt.setBackground(Color.WHITE);
        plt.setOpaque(true);
        plt.setWorldSize(0,0,Math.ceil(frame.kEnd2),0.8);
        plt.setXTicks(10);
        plt.setYTicks(8);
        plt.setXLabel("k/kp");
        plt.setYLabel("ω/ωp");
        plt.setTitle("Dispersion relation plot");
        plt.setSubtitle("");
        
        PlotDialog pltDlg = new PlotDialog(plt,frame,true);
        ArrayList<Double> YPlotA = new ArrayList<>();
        ArrayList<Double> YPlotB = new ArrayList<>();
        ArrayList<Double> YPlotC = new ArrayList<>();
        ArrayList<Double> YPlotD = new ArrayList<>();
        ArrayList<Double> X = new ArrayList<>();

        double wP =  frame.omegaP;
        double kP =  frame.omegaP/c;
       
        Equation eq = new Equation();
        Solver solver = new Solver();
        
        for (double k =frame.kStart1; k<=frame.kEnd1; k+=frame.kStep1) {
            X.add(k); 
            double omegaA = omega(k*kP, frame.omegaP, frame.epsilonHigh);
            double omegaB = omega(k*kP, frame.omegaP, frame.epsilonLow);
            
            eq.setValues(wP, k*kP, frame.epsilonLow, frame.epsilonHigh, frame.thickness);
            double omegaC = Complex.abs(solver.solve(omegaA,omegaB,eq));


            eq.setValues(wP, k*kP, frame.epsilonHigh, frame.epsilonLow, frame.thickness);
            double omegaD = Complex.abs(solver.solve(omegaA,omegaB,eq));

            YPlotA.add(omegaA/wP);
            YPlotB.add(omegaB/wP); 
            YPlotC.add(omegaC/wP); 
            YPlotD.add(omegaD/wP); 
            if (frame.kStep1 == 0.0)
                break;
        }       
        for (double k =frame.kStart2; k<=frame.kEnd2; k+=frame.kStep2) {
            X.add(k); 
            double omegaA = omega(k*kP, frame.omegaP, frame.epsilonHigh);
            double omegaB = omega(k*kP, frame.omegaP, frame.epsilonLow);

            eq.setValues(wP, k*kP, frame.epsilonLow, frame.epsilonHigh, frame.thickness);
            double omegaC = Complex.abs(solver.solve(1.01*omegaA,0.9*omegaB,eq));

            eq.setValues(wP, k*kP, frame.epsilonHigh, frame.epsilonLow, frame.thickness);
            double omegaD = Complex.abs(solver.solve(omegaA,omegaB,eq));

            YPlotA.add(omegaA/wP);
            YPlotB.add(omegaB/wP); 
            YPlotC.add(omegaC/wP); 
            YPlotD.add(omegaD/wP); 
            if (frame.kStep2 == 0.0)
                break;
        }       
        plt.addPlot(X,YPlotA,"A");
        plt.addPlot(X,YPlotB,"B");
        plt.addPlot(X,YPlotC,"C");
        plt.addPlot(X,YPlotD,"D");

        pltDlg.setLocationRelativeTo(null);
        pltDlg.setVisible(true);
    }

}
