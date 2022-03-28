package org.hepforge.alohep;

public class Launcher 
{

    public static void main(final String[] args) 
    {
        final AloHEP alohep = new AloHEP();
        alohep.pack();
        alohep.setLocationRelativeTo(null);
        alohep.setDefaultCloseOperation(3);
        alohep.setVisible(true);
    }
}
