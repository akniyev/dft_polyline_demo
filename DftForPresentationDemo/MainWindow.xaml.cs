using DiscreteFunctionsPlots;
using GraphBuilders;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace DftForPresentationDemo
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        int m;
        double[] vertexes;
        double[] vertexes_y;
        int _N = 16;
        int N
        {
            get
            {
                return _N;
            }
            set
            {
                _N = value;
                lbl_N.Content= _N.ToString();
                if (n > _N) n = _N;
                txt_n.Text = n.ToString();
            }
        }

        int n = 8;

        GraphBuilder2DForm gb = new GraphBuilder2DForm(2);
        DiscreteFunctions.DiscreteFunction2D df;
        Plot2D plot = new Plot2D("DF");
        Plot2D dftPlot = new Plot2D("Fourier");
        Plot2D diffPlot = new Plot2D("Diff");
        Plot2D estPlot = new Plot2D("Estimate");

        public MainWindow()
        {
            InitializeComponent();
            lbl_N.Content = N.ToString();
            gb.Show();
            gb.GraphBuilder.DrawPlot(plot);
            gb.GetGraphBuilder(1).DrawPlot(estPlot);
            gb.GetGraphBuilder(1).DrawPlot(diffPlot);
            gb.GraphBuilder.DrawPlot(dftPlot);
        }

        private void button_Click(object sender, RoutedEventArgs e)
        {
            var tr = new Transformers.FastFourierTransformer();

            var forward = tr.Transform(df.Y.Take(df.Y.Length - 1).ToArray());
            var forward_abs = forward.Select(x => Math.Sqrt(x.x * x.x + x.y * x.y)).ToArray();


            var skipCount = forward.Length / 2 - n;

            for (int i = 0; i < skipCount; i++)
            {
                forward[forward.Length / 2 + i] = 0;
                forward[forward.Length / 2 - i] = 0;
            }

            var backward = tr.InvTransform(forward);
            var backward_list = backward.ToList();
            backward_list.Add(backward[0]);
            backward = backward_list.ToArray();


            var df_dft = ((forward.Length) / 2.0) * new DiscreteFunctions.DiscreteFunctionComplex2D(backward).Re();

            df_dft.X = df.X;

            dftPlot.DiscreteFunction = df_dft;

            dftPlot.Refresh();

            diffPlot.DiscreteFunction = (plot.DiscreteFunction - dftPlot.DiscreteFunction).Abs();

            diffPlot.Refresh();

            //Drawing estimate for S_n(f,x)
            //double[] ksi = new double[m + 1];
            //for (int i = 0; i < m + 1; i++)
            //{
            //    ksi[i] = vertexes[i] - Math.PI;
            //}

            //double[] alpha = new double[m];
            //double[] beta = new double[m];

            //double[] f = new double[m + 1];

            //for (int i = 0; i < m + 1; i++)
            //{
            //    f[i] = vertexes_y[i];
            //}

            //for (int i = 0; i < m; i++)
            //{
            //    alpha[i] = (f[i + 1] - f[i]) / (ksi[i + 1] - ksi[i]);
            //    beta[i] = f[i] - (ksi[i] * (f[i + 1] - f[i]) / (ksi[i + 1] - ksi[i]));
            //}

            //int dN = 10 * N;

            //double[] est = new double[dN];
            //double[] _x = new double[dN];

            //for (int i = 0; i < dN; i++)
            //{
            //    _x[i] = -Math.PI + 2 * Math.PI * i / (dN - 1);
            //}

            //for (int i = 0; i < dN; i++)
            //{
            //    double est_i = 0;

            //    for (int k = n + 1; k < 1000; k++)
            //    {
            //        double inner_sum = 0;
            //        for (int j = 0; j < m; j++)
            //        {
            //            //inner_sum += alpha[j] * Math.Sin(k * ((ksi[j + 1] + ksi[j]) / 2.0 - _x[i])) * Math.Sin(k * ((ksi[j + 1] + ksi[j]) / 2.0));
            //            inner_sum += alpha[j] * (Math.Cos((ksi[j + 1] - _x[i]) * k) - Math.Cos((ksi[j] - _x[i]) * k));
            //        }
            //        inner_sum *= 1.0 / (k * k);

            //        est_i += inner_sum;
            //    }
            //    est_i *= 1.0 / Math.PI;
            //    est[i] = est_i;
            //}

            //double[] _x_shifted = new double[dN];

            //for (int i = 0; i < dN; i++)
            //{
            //    _x_shifted[i] = _x[i] + Math.PI;
            //}
            //var dfEst = new DiscreteFunctions.DiscreteFunction2D(_x_shifted, est);

            //estPlot.DiscreteFunction = dfEst.Abs();

            //estPlot.Refresh();

            //
        }

        private void btn_N_plus_Click(object sender, RoutedEventArgs e)
        {
            if (N < 64000) N *= 2;
        }

        private void btn_N_minus_Click(object sender, RoutedEventArgs e)
        {
            if (N > 4) N /= 2;
        }

        private void textBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            
        }

        private void txt_n_KeyDown(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Enter)
            {
                var b = "";

                for (int i = 0; i < txt_n.Text.Length; i++)
                {
                    if (txt_n.Text[i] <= '9' && txt_n.Text[i] >= '0')
                    {
                        b += txt_n.Text[i];
                    }
                }

                txt_n.Text = b;

                var temp_n = b != "" ? int.Parse(b) : 0;

                if (temp_n < 2) temp_n = 2;

                if (temp_n > N) temp_n = N;

                n = temp_n;
            }
        }

        private void txt_n_LostFocus(object sender, RoutedEventArgs e)
        {
            var b = "";

            for (int i = 0; i < txt_n.Text.Length; i++)
            {
                if (txt_n.Text[i] <= '9' && txt_n.Text[i] >= '0')
                {
                    b += txt_n.Text[i];
                }
            }

            txt_n.Text = b;

            var temp_n = b != "" ? int.Parse(b) : 0;

            if (temp_n < 2) temp_n = 2;

            if (temp_n > N) temp_n = N;

            n = temp_n;
        }

        double[] GenerateGrid(int count)
        {
            var result = new double[count + 1];

            for (int i = 0; i <= count; i++)
            {
                result[i] = 2 * Math.PI * i / N;
            }

            return result;
        }

        private void btn_draw_abs_Click(object sender, RoutedEventArgs e)
        {
            var x = GenerateGrid(N);
            var y = new double[x.Length];

            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Math.Abs(x[i] - Math.PI);
            }

            df = new DiscreteFunctions.DiscreteFunction2D(x, y);

            plot.DiscreteFunction = df;

            plot.Refresh();
        }

        private void btn_draw_sign_Click(object sender, RoutedEventArgs e)
        {
            var x = GenerateGrid(N);
            var y = new double[x.Length];

            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Math.Sign(x[i] - Math.PI);
            }

            df = new DiscreteFunctions.DiscreteFunction2D(x, y);

            plot.DiscreteFunction = df;

            plot.Refresh();
        }

        private void btn_draw_chain_Click(object sender, RoutedEventArgs e)
        {
            m = int.Parse(txt_m.Text);
            var min = int.Parse(txt_min.Text);
            var max = int.Parse(txt_max.Text);

            var x = GenerateGrid(N);
            var y = new double[x.Length];

            vertexes = new double[m + 1];
            var delta = Math.PI * 2.0 / m;
            var r = new Random();
            for (int i = 1; i < m; i++)
            {
                vertexes[i] = delta / 2 + (i - 1) * delta + r.NextDouble() * delta;
            }

            vertexes[0] = 0;
            vertexes[m] = 2 * Math.PI;

            vertexes_y = new double[m + 1];

            for (int i = 0; i < m; i++)
            {
                vertexes_y[i] = r.NextDouble() * (max - min) + min;
            }

            vertexes_y[m] = vertexes_y[0];



            for (int i = 0; i < x.Length - 1; i++)
            {
                int vertexId = 0;
                for (int j = 0; j < m; j++)
                {
                    if (vertexes[j] <= x[i] && vertexes[j + 1] >= x[i])
                    {
                        vertexId = j;
                        break;
                    }
                }

                double a = vertexes[vertexId];
                double b = vertexes[vertexId + 1];

                double c = x[i];

                double f_a = vertexes_y[vertexId];
                double f_b = vertexes_y[vertexId + 1];

                y[i] = f_b * (c - a) / (b - a) + f_a * (b - c) / (b - a);
            }

            y[y.Length - 1] = y[0];

            df = new DiscreteFunctions.DiscreteFunction2D(x, y);

            plot.DiscreteFunction = df;

            plot.Refresh();
        }
    }
}
