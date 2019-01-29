package omegaPlots;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.font.FontRenderContext;
import java.awt.font.TextLayout;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import javax.imageio.ImageIO;
import javax.swing.JComponent;

public class Plot extends JComponent {
    String XLabel;
    String YLabel;
    String title;
    String subtitle;
    double [] WindowLimits;
    double [] WorldLimits;
    double [] BoxLimits;
    double [] BoxInsets;
    ArrayList [] PlotArray;
    int XTicks;
    int YTicks;
    Color [] lineColors;
    class Line {
        double x1;
        double y1;
        double x2;
        double y2;
        Color color;
        String name;
        public Line(double x1, double y1, double x2, double y2,Color col,String n) {
            this.x1 = x1;
            this.y1 = y1;
            this.x2 = x2;
            this.y2 = y2;
            color=col;
            name = n;
        }
    }
    
    public Plot() {
        title = "";
        subtitle = "";
        PlotArray=null;
        WindowLimits=null;
        WorldLimits=null;
        BoxLimits=null;
        BoxInsets=new double[] {70,40,120,70}; // {70,40,30,70};
        XTicks = 5;
        YTicks =5;
        float colors[][] = new float[][]{
            {0F, 0.4470F, 0.7410F},
            {0.8500F, 0.3250F, 0.0980F},
            {0.9290F, 0.6940F, 0.1250F},
            {0.4940F, 0.1840F, 0.5560F},
            {0.4660F, 0.6740F, 0.1880F},
            {0.3010F, 0.7450F, 0.9330F},
            {0.6350F, 0.0780F, 0.1840F},
            {0F, 0.4470F, 0.7410F},
            {0.8500F, 0.3250F, 0.0980F},
            {0.9290F, 0.6940F, 0.1250F},
            {0.4940F, 0.1840F, 0.5560F},
            {0.4660F, 0.6740F, 0.1880F},
            {0.3010F, 0.7450F, 0.9330F},
            {0.6350F, 0.0780F, 0.1840F},
            {0F, 0.4470F, 0.7410F},
            {0.8500F, 0.3250F, 0.0980F},
            {0.9290F, 0.6940F, 0.1250F},
            {0.4940F, 0.1840F, 0.5560F},
            {0.4660F, 0.6740F, 0.1880F},
            {0.3010F, 0.7450F, 0.9330F},
            {0.6350F, 0.0780F, 0.1840F},
            {0F, 0.4470F, 0.7410F},
            {0.8500F, 0.3250F, 0.0980F},
            {0.9290F, 0.6940F, 0.1250F},
            {0.4940F, 0.1840F, 0.5560F},
            {0.4660F, 0.6740F, 0.1880F},
            {0.3010F, 0.7450F, 0.9330F},
            {0.6350F, 0.0780F, 0.1840F},
            {0F, 0.4470F, 0.7410F},
            {0.8500F, 0.3250F, 0.0980F},
            {0.9290F, 0.6940F, 0.1250F},
            {0.4940F, 0.1840F, 0.5560F},
            {0.4660F, 0.6740F, 0.1880F},
            {0.3010F, 0.7450F, 0.9330F},
            {0.6350F, 0.0780F, 0.1840F},
            {0F, 0.4470F, 0.7410F},
            {0.8500F, 0.3250F, 0.0980F},
            {0.9290F, 0.6940F, 0.1250F},
            {0.4940F, 0.1840F, 0.5560F},
            {0.4660F, 0.6740F, 0.1880F},
            {0.3010F, 0.7450F, 0.9330F},
            {0.6350F, 0.0780F, 0.1840F},
            {0F, 0.4470F, 0.7410F},
            {0.8500F, 0.3250F, 0.0980F},
            {0.9290F, 0.6940F, 0.1250F},
            {0.4940F, 0.1840F, 0.5560F},
            {0.4660F, 0.6740F, 0.1880F},
            {0.3010F, 0.7450F, 0.9330F},
            {0.6350F, 0.0780F, 0.1840F},
            {0F, 0.4470F, 0.7410F},
            {0.8500F, 0.3250F, 0.0980F},
            {0.9290F, 0.6940F, 0.1250F},
            {0.4940F, 0.1840F, 0.5560F},
            {0.4660F, 0.6740F, 0.1880F},
            {0.3010F, 0.7450F, 0.9330F},
            {0.6350F, 0.0780F, 0.1840F},
            {0F, 0.4470F, 0.7410F},
            {0.8500F, 0.3250F, 0.0980F},
            {0.9290F, 0.6940F, 0.1250F},
            {0.4940F, 0.1840F, 0.5560F},
            {0.4660F, 0.6740F, 0.1880F},
            {0.3010F, 0.7450F, 0.9330F},
            {0.6350F, 0.0780F, 0.1840F},
            {0F, 0.4470F, 0.7410F}
        };
        lineColors = new Color[colors.length];
        for(int i=0; i<colors.length; i++) {
            lineColors[i] = new Color(colors[i][0],colors[i][1],colors[i][2]);
        }
    }

    public void setTitle(String t) {
        title = t;
    }            

    public void setSubtitle(String t) {
        subtitle = t;
    }            
    
    public void setBoxSize(double x1,double y1,double x2,double y2) {
        BoxLimits=new double[4];
        BoxLimits[0]=x1;
        BoxLimits[1]=y1;
        BoxLimits[2]=x2;
        BoxLimits[3]=y2;
    }

    public void setWindowSize(double x1,double y1,double x2,double y2) {
        WindowLimits=new double[4];
        WindowLimits[0]=x1;
        WindowLimits[1]=y1;
        WindowLimits[2]=x2;
        WindowLimits[3]=y2;
        setBoxSize(x1+BoxInsets[0],y2-BoxInsets[3],x2-BoxInsets[2],y1+BoxInsets[1]);
        setSize((int)(x2-x1),(int)(y2-y1));
    }

    public void setWorldSize(double x1,double y1,double x2,double y2) {
        WorldLimits=new double[4];
        WorldLimits[0]=x1;
        WorldLimits[1]=y1;
        WorldLimits[2]=x2;
        WorldLimits[3]=y2;
    }

    public void setXLabel(String l) {
        XLabel=l;
    }

    public void setYLabel(String l) {
        YLabel=l;
    }

    public double [] transform_coords(double xw,double yw,double [] winLim,double [] wrlLim) {
        double sx,sy;
        double  x1=winLim[0],
                y1=winLim[1],
                x2=winLim[2],
                y2=winLim[3];
        double MIN_X=wrlLim[0];
        double MIN_Y=wrlLim[1];
        double MAX_X=wrlLim[2];
        double MAX_Y=wrlLim[3];
        sx=(double)(x2-x1)/(double)(MAX_X-MIN_X);
        sy=(double)(y2-y1)/(double)(MAX_Y-MIN_Y);
        double xv=(x1+(double)(xw-MIN_X)*sx);
        double yv=(y1+(double)(yw-MIN_Y)*sy);
        return new double[]{xv,yv};
    }
    
    private void drawBox(double [] win,Graphics2D g2) {
        g2.draw(new Line2D.Double(win[0],win[1],win[2],win[1]));
        g2.draw(new Line2D.Double(win[2],win[1],win[2],win[3]));
        g2.draw(new Line2D.Double(win[2],win[3],win[0],win[3]));
        g2.draw(new Line2D.Double(win[0],win[3],win[0],win[1]));
    }

    
    boolean pointInside(Point2D p)
    {
        return (p.getX()>=WorldLimits[0] && p.getX()<=WorldLimits[2] &&
            p.getY()>=WorldLimits[1] && p.getY()<=WorldLimits[3]);
    }

    void clipLine(Line2D line)
    {
      double dx, dy; //, tE, tL;
      Point2D a = line.getP1();
      Point2D b = line.getP2();
      
      dx = b.getX() - a.getX();
      dy = b.getY() - a.getY();

      if (Math.abs(dx)<1e-5 && Math.abs(dy)<1e-5 && pointInside(a)) return;

      Point2D t = new Point2D.Double(0,1);
      
      if (clipT(WorldLimits[0] - a.getX(),  dx, t) &&
          clipT(a.getX() - WorldLimits[2], -dx, t) &&
          clipT(WorldLimits[1] - a.getY(),  dy, t) &&
          clipT(a.getY() - WorldLimits[3], -dy, t))

        {
          if (t.getY() < 1) {
            double bx = a.getX() + t.getY() * dx;
            double by = a.getY() + t.getY() * dy;
            b.setLocation(bx, by);
            line.setLine(a, b);
          }

          if (t.getX() > 0) {
            double ax = a.getX() + t.getX() * dx;
            double ay = a.getY() + t.getX() * dy;
            a.setLocation(ax,ay);
            line.setLine(a, b);
          }

        }
    }

    boolean clipT(double num, double denom, Point2D T) {
      double t;

      if (Math.abs(denom)<1e-5) return (num <= 0.0);

      t = num / denom;

      if (denom > 0) {

        if (t > T.getY()) return false;
        if (t > T.getX()) T.setLocation(t, T.getY());

      } else {

        if (t < T.getX()) return false;
        if (t < T.getY()) T.setLocation(T.getX(), t);

      }

      return true;

    }
    
    
    public void drawWorldLine(double x1,double y1,double x2,double y2,Graphics2D g2) {
        Line2D l = new Line2D.Double(x1,y1,x2,y2);
        if(!pointInside(l.getP1()) && !pointInside(l.getP2()))
            return;
        clipLine(l);
        double [] a=transform_coords(l.getX1(),l.getY1(),BoxLimits,WorldLimits);
        double [] b=transform_coords(l.getX2(),l.getY2(),BoxLimits,WorldLimits);
        Line2D line = new Line2D.Double(a[0],a[1],b[0],b[1]);
        g2.draw(line);
    }

    void drawWorldLine(Line l,Graphics2D g2) {
        drawWorldLine(l.x1,l.y1,l.x2,l.y2,g2);
    }

    public void drawPlots(Graphics2D g2) {
        for(int i=0; i<PlotArray.length; i++) {
            ArrayList pa=PlotArray[i];
            for(int j=0; j<pa.size(); j++) {
                Line [] l=(Line [] )pa.get(j);
                for (int k=0; k<l.length; k++) {
                    g2.setColor(l[k].color);
                    drawWorldLine(l[k],g2);
                }
            }
        }
    }
    
    public void setXTicks(int nticks) {
        XTicks = nticks;
    }

    public void setYTicks(int nticks) {
        YTicks = nticks;
    }
    
    public void draw(Graphics g) {
        Graphics2D g2 = (Graphics2D)g;

        g2.setColor(getBackground());
        g2.fillRect(0, 0, getWidth(), getHeight());
        drawGrid(g2);
        g2.setColor(Color.black);
        g2.setStroke(new BasicStroke());
        drawBox(BoxLimits,g2);
        double wb=BoxLimits[2]-BoxLimits[0];
        double hb=BoxLimits[1]-BoxLimits[3];
        drawStringCentered(1,YLabel,20,hb/2+BoxInsets[1],g2);
        drawStringCentered(0,XLabel,wb/2+BoxInsets[0],BoxInsets[1]+hb+BoxInsets[3]-10,g2);        
        if(!title.isEmpty())
            drawStringCentered(0,title,getWidth()/2,18,new Font("Arial", Font.BOLD, 16),g2);
        if(!subtitle.isEmpty())
            drawStringCentered(0,subtitle,getWidth()/2,34,new Font("Arial", Font.BOLD, 14),g2);
        if (PlotArray!=null) {
            g2.setColor(Color.black);
            g2.setStroke(new BasicStroke());
            drawPlots(g2);
            drawLegend(g2);
        }        
    }
    
    @Override
    public void paintComponent(Graphics g) {
        if (BoxLimits==null) return;
        draw(g);
    }

    protected float [] getStringSize(String s,Font font, Graphics2D g2) {
        FontRenderContext frc = g2.getFontRenderContext();
        Rectangle2D bounds = font.getStringBounds(s, frc);
        float wid = (float)bounds.getWidth();
        float hei = (float)bounds.getWidth();
        return new float[]{wid,hei};
    }

    protected float [] getStringSize(String s, Graphics2D g2) {
        Font font = new Font("Arial", Font.PLAIN, 12);
        return getStringSize(s, font, g2);
    }

    protected void drawString(String s,double x, double y, double theta,Font font,Graphics2D g2) {
        g2.translate(x, y);
        g2.rotate(theta);
        TextLayout textLayout = new TextLayout(s, font,g2.getFontRenderContext());
        textLayout.draw(g2, 0, 0);
        g2.rotate(-theta);
        g2.translate(-x, -y);
    }
    
    protected void drawStringCentered(int dir,String s,double x, double y,Font font,Graphics2D g2) {
      float [] sz=getStringSize(s, font, g2);
      float wid=sz[0];

      if (dir==0) drawString(s, x - wid / 2, y,0,font,g2);
      if (dir==1) drawString(s, x, y+wid/2,-Math.PI/2,font,g2);
    }   

    protected void drawStringCentered(int dir,String s,double x, double y,Graphics2D g2) {
        Font font = new Font("Arial", Font.PLAIN, 12);
        drawStringCentered(dir, s,x, y, font, g2);
    }
    
    private void drawGrid(Graphics2D g2) {
        double firstX = BoxLimits[0];
        double firstY = BoxLimits[1];
        double lastX = BoxLimits[2];
        double lastY = BoxLimits[3];
        g2.setColor(Color.black);

        Stroke stroke1 = new BasicStroke(1,
            BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0,
            new float[] { 2, 3 }, 0);
        Stroke stroke2 = new BasicStroke(1,BasicStroke.CAP_BUTT,BasicStroke.JOIN_BEVEL);

        Line2D line;
        double TickHeight=5;
        drawStringCentered(0, String.format("%5.0g",WorldLimits[0]), firstX, firstY+5*TickHeight, g2);
        drawStringCentered(0, String.format("%5.2g",WorldLimits[2]), lastX, firstY+5*TickHeight, g2);
        String f=String.format("%5.1g",WorldLimits[1]);
        float [] sz=getStringSize(f,g2);
        drawStringCentered(0, f, firstX-TickHeight-sz[0]/2, firstY+4, g2);
        f=String.format("%5.1g",WorldLimits[3]);
        sz=getStringSize(f,g2);
        drawStringCentered(0, f, firstX-TickHeight-sz[0]/2, lastY+4, g2);

        double gridSpaceX=(lastX-firstX)/XTicks;
        double labelSpaceX=(WorldLimits[2]-WorldLimits[0])/XTicks;

        double x = firstX+gridSpaceX;
        double labelx = WorldLimits[0]+labelSpaceX;
        while (x < lastX) {
            drawStringCentered(0, String.format("%5.2g",labelx), x, firstY+5*TickHeight, g2);
            line = new Line2D.Double(x, firstY, x, lastY);
            g2.setStroke(stroke1);
            g2.draw(line);
            g2.setStroke(stroke2);
            line = new Line2D.Double(x, firstY, x, firstY-TickHeight);
            g2.draw(line);
            line = new Line2D.Double(x, lastY, x, lastY+TickHeight);
            g2.draw(line);
            x += gridSpaceX;
            labelx += labelSpaceX;
        }

        double gridSpaceY=(-lastY+firstY)/YTicks;
        double labelSpaceY=(WorldLimits[3]-WorldLimits[1])/YTicks;

        double y = firstY-gridSpaceY;
        double labely = WorldLimits[1]+labelSpaceY;
        while (y > lastY) {
            f=String.format("%5.1g",labely);
            sz=getStringSize(f,g2);
            drawStringCentered(0, f, firstX-TickHeight-sz[0]/2, y+4, g2);
            line = new Line2D.Double(firstX, y, lastX, y);
            g2.setStroke(stroke1);
            g2.draw(line);
            g2.setStroke(stroke2);
            line = new Line2D.Double(firstX, y, firstX+TickHeight, y);
            g2.draw(line);
            line = new Line2D.Double(lastX, y, lastX-TickHeight, y);
            g2.draw(line);
            y -= gridSpaceY;
            labely += labelSpaceY;
        }
    }
    
    public int addPlot(){
        if(PlotArray==null) {
            PlotArray=new ArrayList[1];
            PlotArray[0]=new ArrayList();
        }else{
            ArrayList [] tmp=PlotArray;
            PlotArray=new ArrayList[PlotArray.length+1];
            int i;
            for(i=0; i<tmp.length; i++){
                PlotArray[i]=tmp[i];
            }
            PlotArray[i]=new ArrayList();
        }
        return (PlotArray.length-1);
    }
   
    public int addPlot(ArrayList<Double> X,ArrayList<Double> Y,String name) {
        int j=addPlot();
        Line [] l=new Line[X.size()-1];
        for (int i=0; i<X.size()-1; i++) {
            l[i]=new Line(X.get(i),Y.get(i),X.get(i+1),Y.get(i+1),lineColors[j%lineColors.length],name);
        }
        PlotArray[j].add(l);
        return j;
    }
    
    public void drawLegend(Graphics2D g2) {
        double [] lbox=new double[4];
        lbox[0]=BoxLimits[2] + 20; // getWidth()*9.3/16;
        lbox[1]=BoxInsets[1]; // getHeight()*3/16;
        lbox[2]=lbox[0] + 80; //260; //getWidth()*14/16;
        lbox[3]=lbox[1]+20*PlotArray.length;
        g2.setBackground(Color.white);
        g2.clearRect((int)lbox[0],(int)lbox[1],(int)(lbox[2]-lbox[0]),(int)(lbox[3]-lbox[1]));
        g2.setColor(Color.black);
        g2.setStroke(new BasicStroke());
        drawBox(lbox,g2);
        Font font = new Font("Arial", Font.PLAIN, 12);
        for (int i=0; i<PlotArray.length; i++) {
            Line [] l=(Line [] )PlotArray[i].get(0);
            double y = lbox[1]+10+20*i;
            g2.setColor(l[0].color);
            g2.draw(new Line2D.Double(lbox[0]+5,y,lbox[0]+35,y));
            g2.setColor(Color.black);
            drawString(l[0].name,lbox[0]+45, y+5, 0, font,g2);
        }
    }
    
    public boolean save(String filename) {
        try {
            BufferedImage img = new BufferedImage(getWidth(), getHeight(),
                    BufferedImage.TYPE_INT_RGB);
            Graphics g = img.getGraphics();
            draw(g);
            ImageIO.write(img, "png", new File(filename));
            return true;
        } catch (IOException e) {
            return false;
        }        
    }
    
    public boolean saveText(String filename) {
        BufferedWriter  fidw;
        try {
            fidw=new BufferedWriter(new FileWriter(filename));
            ArrayList lines = new ArrayList();
            String header = "";
            for(int i=0; i<PlotArray.length; i++) {
                Line [] l=(Line [] )PlotArray[i].get(0);
                String name = l[0].name.replace(",", "");
                if(i>0)
                    header+=",";
                header+=String.format("%s",name);
            }
            fidw.write(this.XLabel + ",");
            fidw.write(header);
            fidw.newLine();
            
            int k=0;
            boolean done = false;            
            while(!done) {
                for(int i=0; i<PlotArray.length; i++) {
                    Line [] l=(Line [] )PlotArray[i].get(0);
                    if(i==0) {
                        if(k < l.length)
                            fidw.write(String.format("%5.6f",l[k].x1));
                        else
                            fidw.write(String.format("%5.6f",l[l.length-1].x2));
                    }
                    fidw.write(",");
                    if(k < l.length)
                        fidw.write(String.format("%5.6f",l[k].y1));
                    else {
                        fidw.write(String.format("%5.6f",l[l.length-1].y2));
                        done = true;
                    }
                }
                fidw.newLine();
                k++;
            }
            fidw.close();
            return true;
        }
        catch (IOException ie){
            return false;
        }        
    }
}
