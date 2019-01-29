package omegaPlots;

import java.io.File;
import javax.swing.filechooser.FileFilter;

class myFileFilter extends FileFilter {
    String[] extensions;
    String description;

    public myFileFilter(String ext,String descr) {
        extensions=new String[1];
        extensions[0]=ext.toLowerCase();
        description=(descr==null? ext + " files":descr);
    }

    @Override
    public boolean accept(File f) {
        if(f.isDirectory()) { return true; }

        String name = f.getName().toLowerCase();
        for (int i = extensions.length - 1; i >= 0; i--) {
            if (name.endsWith("."+extensions[i])) {
                return true;
            }
        }
        return false;
    }

    @Override
    public String getDescription() { 
        return description; 
    }
}    
