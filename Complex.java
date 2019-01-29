package omegaPlots;

public class Complex {
   public double  re;
   public double  im;
   
   public Complex(double r,double i) {
        re = r;
        im = i;
   }
   public Complex(double r) {
       re = r;
       im = 0.0;
   }
   public static double real(Complex a)
   {
       return a.re;
   }
   public static double imag(Complex a)
   {
       return a.im;
   }
   public static Complex I()
   {
       return new Complex(0.0,1.0);
   }
   public static Complex add(Complex a,Complex b) {
        double re = a.re + b.re;
        double im= a.im + b.im;
        return new Complex(re,im);
   }

   public static Complex add(Complex a,double b) {
       return add(a,new Complex(b));
   }
   
   public static Complex sub(Complex a,Complex b) {
        double re = a.re - b.re;
        double im= a.im - b.im;
        return new Complex(re,im);
   }
   
   public static Complex mul(Complex a,Complex b) {
        double re = a.re*b.re - a.im*b.im;
        double im= a.re*b.im + a.im*b.re;
        return new Complex(re,im);
   }

   public static Complex mul(Complex a,double b) {
        double re = a.re*b;
        double im= a.im*b;
        return new Complex(re,im);
   }
   public static Complex div(Complex a,Complex b) {
       double bb = Math.pow(Complex.abs(b),2);
       double re = (a.re*b.re + a.im*b.im)/bb;
       double im = (a.im*b.re - a.re*b.im)/bb;
       return new Complex(re,im);
   }
   public static double abs(Complex a) {      
        return Math.hypot(a.re,a.im);          
   }
   
   public static Complex sqrt(Complex a) {       
        double re = Math.sqrt((Complex.abs(a)+a.re)/2.0);
        double im = Math.sqrt((Complex.abs(a)-a.re)/2.0);
        double sign;
        sign =(a.im>=0)?  1.0 : -1.0;
        return new Complex(re,sign*im);
   }

   public static Complex sqrt(double re) {       
        return sqrt(new Complex(re));
   }

   public static Complex log(Complex a) {
       double re = Math.log(Complex.abs(a));
       double im = Math.atan2(a.im,a.re);
       return new Complex(re,im);
   }
   
   public static Complex asin(Complex z) {
       Complex iz = new Complex(-z.im,z.re);
       Complex z2 = mul(z, z);
       Complex sq1_z2 = sqrt(sub(new Complex(1),z2));
       Complex arg = Complex.add(iz,sq1_z2);
       return Complex.div(Complex.log(arg),Complex.I());
   }

   public static Complex sin(Complex a) {
       if(a.im==0)
           return new Complex(Math.sin(a.re));
       double re = Math.sin(a.re)*Math.cosh(a.im);
       double im = Math.cos(a.re)*Math.sinh(a.im);
       return new Complex(re,im);
   }

   public static Complex cos(Complex a) {
       if(a.im==0)
           return new Complex(Math.cos(a.re));
       double re = Math.cos(a.re)*Math.cosh(a.im);
       double im = -Math.sin(a.re)*Math.sinh(a.im);
       return new Complex(re,im);
   }
   public static Complex tan(Complex z) {
       return div(sin(z),cos(z));
   }
   public static Complex exp(Complex z) {
       double ex = Math.exp(z.re);
       double re = ex*Math.cos(z.im);
       double im = ex*Math.sin(z.im);
       return new Complex(re,im);
   }
   public static Complex conj(Complex z) {
       return new Complex(z.re,-z.im);
   }

    @Override
    public String toString() {
        if(Math.abs(im)<1e-10)
            return ""+re;
        else
            return re+((im>0)? "+":"-")+Math.abs(im)+"i";
    }

    public static Complex tanh(Complex z) {
       if(z.im==0)
            return new Complex(Math.tanh(z.re));
       else
           return div(tan(mul(I(),z)),I());
    }

}
