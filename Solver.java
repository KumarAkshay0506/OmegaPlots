package omegaPlots;

public class Solver {
    double MAXERROR;
    double MAXITERATIONS;
    
    public Solver() {
        MAXERROR = 1e-6;
        MAXITERATIONS = 200;
    }
    
    public Solver(double maxerror, int maxiterations) {
        MAXERROR = maxerror;
        MAXITERATIONS = maxiterations;
    }

    Complex solve(double x0,double x1,Function f)
    {
        Complex x, _x0, _x1;
        int iterations;

        _x0 = new Complex(x0);
        _x1 = new Complex(x1);
        Complex fx0 = f.eval(_x0);
        Complex fx1 = f.eval(_x1);
        x = Complex.div(Complex.sub(Complex.mul(_x0,fx1), Complex.mul(_x1,fx0)),Complex.sub(fx1, fx0));   
        iterations=0;
        while(Math.abs(Complex.real(f.eval(x))) > MAXERROR)
        {
            _x0 = _x1;
            _x1 = x;
            fx0 = fx1;
            fx1 = f.eval(_x1);
            x = Complex.div(Complex.sub(Complex.mul(_x0,fx1), Complex.mul(_x1,fx0)),Complex.sub(fx1, fx0)); 
            iterations++;
            if(iterations>=MAXITERATIONS)
                break;
        }
        return x;
    }    
    
}
