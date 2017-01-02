package hacettepe.lgmd.lrohs;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Graphics;

import javax.swing.JFrame;
import javax.swing.JPanel;

public class TestPanel extends JPanel {
	
   	private static final long serialVersionUID = 1L;
   	
   	//read http://stackoverflow.com/questions/15544549/how-does-paintcomponent-work
//   	http://www.oracle.com/technetwork/java/painting-140037.html

	public void paintComponent(Graphics g) {
        super.paintComponent(g);
        
        g.setColor(Color.BLACK);
        g.drawRect(10,10, 50,250);
        
        g.setColor(Color.RED);
        g.drawLine(5, 50, 65, 50);
        
        g.setColor(Color.BLACK);
        g.fillRect(10,15, 50,40);

        g.setColor(Color.BLACK);
        g.fillRect(250,100, 20,20);

//        for (int i = 0; i < 20; i++) {
//            draw(g);
//        }
    }

    public void draw(Graphics g) {
        Color c = new Color((int) (Math.random() * 255), (int) (Math.random() * 255), (int) (Math.random() * 255));
        g.setColor(c);
        g.fillRect((int) (Math.random() * 400), (int) (Math.random() * 300), (int) (Math.random() * 40), (int) (Math.random() * 40));
    }

    public static void main(String[] args) {
        JFrame f = new JFrame();
        f.getContentPane().add(new TestPanel(), BorderLayout.CENTER);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setSize(400, 300);
        f.setVisible(true);
    }
}