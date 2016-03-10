using DiscreteFunctions;
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
using static System.Math;

namespace Polyline
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        GraphBuilder2DForm gb = new GraphBuilder2DForm(2);
        Plot2D plot_polyline = new Plot2D("Polyline");
        Plot2D plot_discrete_polyline = new Plot2D("Discrete Polyline");
        Plot2D plot_fourier = new Plot2D("Fourier");
        Plot2D plot_dft = new Plot2D("DFT");
        Plot2D plot_fourier_diff = new Plot2D("Fourier diff");
        Plot2D plot_calculated_fourier_diff = new Plot2D("Calculated Fourier Diff");
        Plot2D plot_dft_diff = new Plot2D("DFT Diff");

        Plot2D plot_uniform_estimate = new Plot2D("Uniform diff estimate");

        //Polyline params
        double[] x_polyline;
        double[] y_polyline;

        double[] x_discrete_polyline;
        double[] y_discrete_polyline;

        double[] alpha_polyline;
        double[] beta_polyline;

        int m;
        double min_y;
        double max_y;


        public MainWindow()
        {
            InitializeComponent();
            gb.GraphBuilder.DrawPlot(plot_polyline);
            gb.GraphBuilder.DrawPlot(plot_discrete_polyline);
            gb.GraphBuilder.DrawPlot(plot_fourier);
            gb.GraphBuilder.DrawPlot(plot_dft);
            gb.GetGraphBuilder(1).DrawPlot(plot_fourier_diff);
            gb.GetGraphBuilder(1).DrawPlot(plot_calculated_fourier_diff);
            gb.GetGraphBuilder(1).DrawPlot(plot_dft_diff);
            gb.GetGraphBuilder(1).DrawPlot(plot_uniform_estimate);
            gb.Show();
        }

        double a0_2()
        {
            double[] ksi = x_polyline;

            double a0_2 = 0;
            for (int j = 0; j < m; j++)
            {
                a0_2 += alpha_polyline[j] / 2.0 * (ksi[j + 1] * ksi[j + 1] - ksi[j] * ksi[j]) + beta_polyline[j] * (ksi[j + 1] - ksi[j]);
            }
            a0_2 *= 1.0 / (Math.PI * 2);

            return a0_2;
        }

        double a_(int k)
        {
            double[] ksi = x_polyline;
            double a_k = 0;
            for (int j = 0; j < m; j++)
            {
                a_k += alpha_polyline[j] * (Math.Cos(k * ksi[j + 1]) - Math.Cos(k * ksi[j]));
            }

            a_k *= 1.0 / (Math.PI * k * k);

            return a_k;
        }

        double b_(int k)
        {
            double[] ksi = x_polyline;
            double b_k = 0;
            for (int j = 0; j < m; j++)
            {
               b_k += alpha_polyline[j] * (Math.Sin(k * ksi[j + 1]) - Math.Sin(k * ksi[j]));
            }

            b_k *= 1.0 / (Math.PI * k * k);

            return b_k;
        }

        void generatePolyline()
        {
            m = int.Parse(txt_m.Text);
            min_y = double.Parse(txt_min_y.Text);
            max_y = double.Parse(txt_max_y.Text);

            var r = new Random();

            x_polyline = new double[m - 1].Select(x => r.NextDouble() * 2 * Math.PI - Math.PI).Union(new double[] { Math.PI, -Math.PI }).OrderBy(y => y).ToArray();
            y_polyline = x_polyline.Select(x => r.NextDouble() * (max_y - min_y) + min_y).ToArray();
            y_polyline[0] = y_polyline[m];

            plot_polyline.DiscreteFunction = new DiscreteFunction2D(x_polyline, y_polyline);
            plot_polyline.Refresh();

            //Computing alpha and beta
            alpha_polyline = new double[m];
            beta_polyline = new double[m];

            for (int i = 0; i < m; i++)
            {
                double x1, x2, y1, y2;
                x1 = x_polyline[i];
                x2 = x_polyline[i + 1];
                y1 = y_polyline[i];
                y2 = y_polyline[i + 1];

                alpha_polyline[i] = (y2 - y1) / (x2 - x1);
                beta_polyline[i] = y1 - alpha_polyline[i] * x1;
            }
        }

        private void btn_generate_polyline_Click(object sender, RoutedEventArgs e)
        {
            generatePolyline();
        }

        private void btn_generate_S_n_Click(object sender, RoutedEventArgs e)
        {
            //Generating discrete polyline
            int N = int.Parse(txt_DFT_N.Text);
            int currentSegment = 0;
            x_discrete_polyline = new double[N];
            y_discrete_polyline = new double[N];

            for (int i = 0; i < N; i++)
            {
                double x = Math.PI * 2.0 * i / N - Math.PI;
                x_discrete_polyline[i] = x;
                while (x > x_polyline[currentSegment + 1])
                    currentSegment++;
                y_discrete_polyline[i] = alpha_polyline[currentSegment] * x + beta_polyline[currentSegment];
            }

            plot_discrete_polyline.DiscreteFunction = new DiscreteFunction2D(x_discrete_polyline, y_discrete_polyline);
            plot_discrete_polyline.Refresh();

            //Checking Fourier series
            //Вычисляем частичную сумму ряда Фурье
            int Dens = int.Parse(txt_Dens.Text);
            int n_fourier = int.Parse(txt_n.Text);

            double[] x_fourier = new double[Dens + 1];
            double[] y_fourier = new double[Dens + 1];
            double[] ksi = x_polyline;

            for (int i = 0; i <= Dens; i++)
            {

                double x = 2.0 * Math.PI / Dens * i - Math.PI;
                x_fourier[i] = x;
                double y = 0;

                //Здесь начинается новый способ вычисления (выведенный по формулам)
                //Новый способ тоже работает, как и старый
                double a0 = 0;
                for (int j = 0; j < m; j++)
                {
                    a0 += alpha_polyline[j] / 2.0 * (ksi[j+1] * ksi[j+1] - ksi[j] * ksi[j]) + beta_polyline[j] * (ksi[j+1] - ksi[j]);
                    
                }
                a0 *= 1.0 / (2.0 * Math.PI);
                y += a0;

                double sum_k = 0;
                for (int k = 1; k <= n_fourier; k++)
                {
                    double sum_j = 0;
                    for (int j = 0; j < m; j++)
                    {
                        sum_j += alpha_polyline[j] * (Math.Cos(k * x) * (Math.Cos(k * ksi[j + 1]) - Math.Cos(k * ksi[j])) + Math.Sin(k * x) * (Math.Sin(k * ksi[j + 1]) - Math.Sin(k * ksi[j])));
                    }
                    sum_j *= 1.0 / (k * k);

                    sum_k += sum_j;
                }
                sum_k *= 1.0 / Math.PI;

                y += sum_k;

                //Это старый способ вычисления, вместо него попробуем новый
                //y = a0_2();
                //for (int k = 1; k <= n; k++)
                //{
                //    y += a_(k) * Math.Cos(k * x) + b_(k) * Math.Sin(k * x);
                //}
                y_fourier[i] = y;

            }
            plot_fourier.DiscreteFunction = new DiscreteFunction2D(x_fourier, y_fourier);
            plot_fourier.Refresh();

            //Вычисляем дискретную сумму ряда Фурье на густой сетке
            int n = int.Parse(txt_DFT_n.Text);

            var tr = new Transformers.FastFourierTransformer();
            var forward = tr.Transform(y_discrete_polyline);

            var forward_partial = new alglib.complex[forward.Length];

            for (int i = 0; i < forward_partial.Length; i++)
            {
                var deg = Min(i, N - i);
                forward_partial[i] = deg <= n ? forward[i] : new alglib.complex(0);
            }

            var left = forward_partial.Take((N + 1)/2).ToArray();
            var right = forward_partial.Skip((N + 1) / 2).ToArray();

            var zeros = new List<alglib.complex>();
            var dens = 1;
            while (dens < Dens) dens *= 2;

            var zerosCount = dens - forward.Length;
            while (zeros.Count < zerosCount) zeros.Add(new alglib.complex(0, 0));

            var forwardDense_list = new List<alglib.complex>();
            forwardDense_list.AddRange(left);
            forwardDense_list.AddRange(zeros);
            forwardDense_list.AddRange(right);

            var forwardDense = forwardDense_list.Select(x => x * dens / 2).ToArray();

            
            
            var y_inverse = tr.InvTransform(forwardDense).Select(x => x.x).ToArray();

            double[] x_inverse = new double[y_inverse.Length];

            for (int i = 0; i < y_inverse.Length; i++)
            {
                x_inverse[i] = 2.0 * PI * i / y_inverse.Length - PI;
            }

            plot_dft.DiscreteFunction = new DiscreteFunction2D(x_inverse, y_inverse);
            plot_dft.Refresh();

            //Diff df - dft
            double[] dft_diff = new double[y_inverse.Length];
            currentSegment = 0;

            for (int i = 0; i < dft_diff.Length; i++)
            {
                double x = x_inverse[i];
                while (x > x_polyline[currentSegment + 1])
                    currentSegment++;
                double y = alpha_polyline[currentSegment] * x + beta_polyline[currentSegment];
                dft_diff[i] = y - y_inverse[i];
            }

            plot_dft_diff.DiscreteFunction = new DiscreteFunction2D(x_inverse, dft_diff);
            plot_dft_diff.Refresh();


            //Diff f - fourier
            //Берем разность функции и вычисленного выше ряда Фурье
            double[] diff = new double[Dens + 1];
            currentSegment = 0;

            for (int i = 0; i <= Dens; i++)
            {
                double x = Math.PI * 2.0 * i / Dens - Math.PI;
                while (x > x_polyline[currentSegment + 1])
                    currentSegment++;
                double y = alpha_polyline[currentSegment] * x + beta_polyline[currentSegment];
                diff[i] = y - y_fourier[i];
            }

            plot_fourier_diff.DiscreteFunction = new DiscreteFunction2D(x_fourier, diff);
            plot_fourier_diff.Refresh();

            //Calculated Fourier Diff
            //Вычисляем самостоятельно остаток ряда Фурье (естественно, не до бесконечности, а до 1000 или около того)

            int j_skipped = int.Parse(txt_j_to_skip.Text);

            double[] calculated_diff_x = new double[Dens + 1];
            double[] calculated_diff_y = new double[Dens + 1];
            //for (int i = 0; i <= N; i++)
            //{
            //    double x = Math.PI * 2.0 * i / N - Math.PI;
            //    calculated_diff_x[i] = x;
            //    double sum_k = 0;
            //    for (int k = n + 1; k < 1000; k++)
            //    {
            //        double sum_j = 0;
            //        for (int j = 1; j <= m; j++)
            //        {
            //            int _j = j;
            //            int _j_1 = j - 1;

            //            if (_j == m) _j = 0;

            //            if (j_skipped != -1 && j != j_skipped) continue;
            //            sum_j += (alpha_polyline[_j_1] - alpha_polyline[_j]) * Math.Cos(k * (ksi[_j] - x));
            //        }
            //        //sum_j += (alpha_polyline[m - 1] - alpha_polyline[0]) * Math.Cos(k * (x + Math.PI));
            //        sum_j *= 1.0 / (k * k);

            //        sum_k += sum_j;
            //    }
            //    sum_k *= 1.0 / Math.PI;
            //    calculated_diff_y[i] = sum_k;
            //}

            //Working, but commented
            for (int i = 0; i <= Dens; i++)
            {
                double x = Math.PI * 2.0 * i / Dens - Math.PI;
                calculated_diff_x[i] = x;
                double sum_j = 0;
                for (int j = 1; j <= m; j++)
                {
                    if (j_skipped != -1 && j != j_skipped) continue;
                    double sum_k = 0;
                    for (int k = n_fourier + 1; k <= 1000 + n_fourier + 1; k++)
                    {
                        sum_k += Math.Cos(k * (ksi[j] - x)) / (k * k);
                    }
                    //sum_j += (alpha_polyline[m - 1] - alpha_polyline[0]) * Math.Cos(k * (x + Math.PI));
                    double alpha_j = alpha_polyline[j % m];
                    double alpha_j_1 = alpha_polyline[(j - 1) % m];
                    sum_k *= alpha_j_1 - alpha_j;

                    sum_j += sum_k;
                }
                sum_j *= 1.0 / Math.PI;
                calculated_diff_y[i] = sum_j;
            }

            plot_calculated_fourier_diff.DiscreteFunction = new DiscreteFunction2D(calculated_diff_x, calculated_diff_y);
            plot_calculated_fourier_diff.Refresh();

            //Здесь вычисляем вычисленный остаток для дискретных рядов Фурье
            //int N = int.Parse(txt_DFT_N.Text);
            //int dft_n = int.Parse(txt_DFT_n.Text);

            //double[] calculated_dft_diff_x = new double[Dens + 1];
            //double[] calculated_dft_diff_y = new double[Dens + 1];

            //for (int i = 0; i <= Dens; i++)
            //{
            //    double x = Math.PI * 2.0 * i / Dens - Math.PI;
            //    calculated_dft_diff_x[i] = x;

            //    double R1 = 0;
            //    for (int mu = 1; mu < 1000; mu++)
            //    {
            //        double sum_j = 0;
            //        for (int j = 1; j <= m; j++)
            //        {
            //            int _j = j % m;
            //            sum_j += (alpha_polyline[j - 1] - alpha_polyline[_j]) * Math.Cos(mu * Dens * ksi[j]);
            //        }
            //        sum_j *= 1.0 / (mu * mu);

            //        R1 += sum_j;
            //    }
            //    R1 *= 1.0 / (Math.PI * N * N);

            //    double R2 = 0;

            //    for (int mu = 0; mu < 1000; mu++)
            //    {
            //        double sum_k = 0;
            //        for (int k = 1; k <= dft_n; k++)
            //        {
            //            double sum_j = 0;
            //            for (int j = 1; j<= m; j++)
            //            {
            //                sum_j += (alpha_polyline[j - 1] - alpha_polyline[j % m]) * ((Math.Cos(k * (x - ksi[j]) + mu * N * ksi[j])) / (Math.Pow(k - mu * N, 2) + (Math.Cos(k * (x - ksi[j]) - mu * N * ksi[j])) / (Math.Pow(k + mu * N, 2))));
            //            }
            //            sum_k += sum_j;
            //        }
            //        R2 += sum_k;
            //    }
            //    R2 *= 1.0 / Math.PI;

            //    double R = R1 + R2;

            //    calculated_dft_diff_y[i] = R;
            //}

            //plot_dft_diff.DiscreteFunction = new DiscreteFunction2D(calculated_dft_diff_x, calculated_dft_diff_y);
            //plot_dft_diff.Refresh();

            //Calculating uniform estimate
            double[] calculated_diff_esimate_uniform_x = new double[2];
            double[] calculated_diff_estimate_uniform_y = new double[2];

            calculated_diff_esimate_uniform_x[0] = -Math.PI;
            calculated_diff_esimate_uniform_x[1] = Math.PI;

            double est = 1.0 / Math.PI / n_fourier;
            double est_sum_j = 0;
            for (int j = 1; j < m; j++)
            {
                est_sum_j += Math.Abs(alpha_polyline[j-1] - alpha_polyline[j]);
            }
            est_sum_j += Math.Abs(alpha_polyline[m-1] - alpha_polyline[0]);

            est *= est_sum_j;

            calculated_diff_estimate_uniform_y[0] = est;
            calculated_diff_estimate_uniform_y[1] = est;

            plot_uniform_estimate.DiscreteFunction = new DiscreteFunction2D(calculated_diff_esimate_uniform_x, calculated_diff_estimate_uniform_y);
            plot_uniform_estimate.Refresh();
        }
    }
}
