import java.util.Scanner;

public class Main {

    public static class ComplexNumber {
        private double re;
        private double im;

        public ComplexNumber(double re, double im) {
            this.re = re;
            this.im = im;
        }

        public ComplexNumber(double re) {
            this(re, 0.0);
        }

        public ComplexNumber add(ComplexNumber other) {
            return new ComplexNumber(this.re + other.re, this.im + other.im);
        }

        public ComplexNumber subtract(ComplexNumber other) {
            return new ComplexNumber(this.re - other.re, this.im - other.im);
        }

        public ComplexNumber multiply(ComplexNumber other) {
            double newRe = this.re * other.re - this.im * other.im;
            double newIm = this.re * other.im + this.im * other.re;
            return new ComplexNumber(newRe, newIm);
        }

        public ComplexNumber divide(ComplexNumber other) {
            double denominator = other.re * other.re + other.im * other.im;
            if (Math.abs(denominator) < 1e-12) {
                throw new ArithmeticException("Деление комплексного числа на ноль");
            }
            double newRe = (this.re * other.re + this.im * other.im) / denominator;
            double newIm = (this.im * other.re - this.re * other.im) / denominator;
            return new ComplexNumber(newRe, newIm);
        }

        @Override
        public String toString() {
            if (Math.abs(im) < 1e-9) return String.format("%.2f", re);
            if (Math.abs(re) < 1e-9) return String.format("%.2fi", im);
            return String.format("%.2f %s %.2fi", re, (im < 0 ? "-" : "+"), Math.abs(im));
        }
    }

    public static class Matrix {
        private final int rows;
        private final int cols;
        private final ComplexNumber[][] data;

        public Matrix(int rows, int cols) {
            if (rows <= 0 || cols <= 0) {
                throw new IllegalArgumentException("Размеры матрицы должны быть больше 0");
            }
            this.rows = rows;
            this.cols = cols;
            this.data = new ComplexNumber[rows][cols];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    data[i][j] = new ComplexNumber(0, 0);
                }
            }
        }

        public int getRows() { return rows; }
        public int getCols() { return cols; }

        public void setElement(int r, int c, ComplexNumber val) {
            data[r][c] = val;
        }

        public ComplexNumber getElement(int r, int c) {
            return data[r][c];
        }

        public Matrix add(Matrix other) {
            if (this.rows != other.rows || this.cols != other.cols) {
                throw new IllegalArgumentException("Матрицы должны быть одинакового размера");
            }
            Matrix result = new Matrix(rows, cols);
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    result.setElement(i, j, this.data[i][j].add(other.data[i][j]));
                }
            }
            return result;
        }

        public Matrix multiply(Matrix other) {
            if (this.cols != other.rows) {
                throw new IllegalArgumentException("Число столбцов первой матрицы должно равняться числу строк второй матрицы");
            }
            Matrix result = new Matrix(this.rows, other.cols);
            for (int i = 0; i < this.rows; i++) {
                for (int j = 0; j < other.cols; j++) {
                    ComplexNumber sum = new ComplexNumber(0, 0);
                    for (int k = 0; k < this.cols; k++) {
                        sum = sum.add(this.data[i][k].multiply(other.data[k][j]));
                    }
                    result.setElement(i, j, sum);
                }
            }
            return result;
        }

        public Matrix transpose() {
            Matrix result = new Matrix(cols, rows);
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    result.setElement(j, i, this.data[i][j]);
                }
            }
            return result;
        }

        public ComplexNumber determinant() {
            if (rows != cols) {
                throw new IllegalArgumentException("Определитель есть только у квадратных матриц");
            }
            return calcDeterminant(this);
        }

        private ComplexNumber calcDeterminant(Matrix m) {
            int n = m.rows;
            if (n == 1) return m.data[0][0];
            if (n == 2) {
                return m.data[0][0].multiply(m.data[1][1])
                        .subtract(m.data[0][1].multiply(m.data[1][0]));
            }

            ComplexNumber det = new ComplexNumber(0, 0);
            for (int j = 0; j < n; j++) {
                Matrix subMatrix = getSubMatrix(m, 0, j);
                ComplexNumber sign = ((j % 2 == 0) ? new ComplexNumber(1) : new ComplexNumber(-1));
                ComplexNumber term = m.data[0][j].multiply(sign).multiply(calcDeterminant(subMatrix));
                det = det.add(term);
            }
            return det;
        }

        private Matrix getSubMatrix(Matrix src, int excludingRow, int excludingCol) {
            Matrix sub = new Matrix(src.rows - 1, src.cols - 1);
            int r = -1;
            for (int i = 0; i < src.rows; i++) {
                if (i == excludingRow) continue;
                r++;
                int c = -1;
                for (int j = 0; j < src.cols; j++) {
                    if (j == excludingCol) continue;
                    c++;
                    sub.setElement(r, c, src.data[i][j]);
                }
            }
            return sub;
        }

        public Matrix divide(Matrix other) {
            Matrix inverse = other.inverse();
            return this.multiply(inverse);
        }

        public Matrix inverse() {
            ComplexNumber det = this.determinant();
            if (Math.abs(det.re) < 1e-12 && Math.abs(det.im) < 1e-12) {
                throw new ArithmeticException("Обратная матрица не существует, так как определитель равен 0");
            }

            Matrix adjugate = new Matrix(rows, cols);
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    Matrix sub = getSubMatrix(this, i, j);
                    ComplexNumber minorDet = sub.calcDeterminant(sub);
                    ComplexNumber sign = (((i + j) % 2 == 0) ? new ComplexNumber(1) : new ComplexNumber(-1));
                    adjugate.setElement(i, j, minorDet.multiply(sign));
                }
            }

            Matrix adjTransposed = adjugate.transpose();
            Matrix inv = new Matrix(rows, cols);
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    inv.setElement(i, j, adjTransposed.getElement(i, j).divide(det));
                }
            }
            return inv;
        }

        public void print() {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    System.out.print("[" + data[i][j] + "]\t");
                }
                System.out.println();
            }
        }
    }

    private static Matrix readMatrix(Scanner sc, String name) {
        System.out.print("Введите число строк матрицы " + name + ": ");
        int n = sc.nextInt();
        System.out.print("Введите число столбцов матрицы " + name + ": ");
        int m = sc.nextInt();

        Matrix matrix = new Matrix(n, m);
        System.out.println("Заполнение матрицы " + name + " (вводите действительную и мнимую часть):");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                System.out.printf("[%d][%d] Действительная часть: ", i + 1, j + 1);
                double re = sc.nextDouble();
                System.out.printf("[%d][%d] Мнимая часть: ", i + 1, j + 1);
                double im = sc.nextDouble();
                matrix.setElement(i, j, new ComplexNumber(re, im));
            }
        }
        return matrix;
    }

    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);

        while (true) {
            try {
                System.out.println("1. Сложить матрицы");
                System.out.println("2. Перемножить матрицы");
                System.out.println("3. Разделить матрицы");
                System.out.println("4. Транспонировать матрицу");
                System.out.println("5. Найти определитель матрицы");
                System.out.println("0. Выход");
                System.out.print("Выберите действие: ");

                int choice = sc.nextInt();
                if (choice == 0) break;

                switch (choice) {
                    case 1: {
                        Matrix m1 = readMatrix(sc, "A");
                        Matrix m2 = readMatrix(sc, "B");
                        System.out.println("Результат сложения:");
                        m1.add(m2).print();
                        break;
                    }
                    case 2: {
                        Matrix m1 = readMatrix(sc, "A");
                        Matrix m2 = readMatrix(sc, "B");
                        System.out.println("Результат умножения:");
                        m1.multiply(m2).print();
                        break;
                    }
                    case 3: {
                        Matrix m1 = readMatrix(sc, "A");
                        Matrix m2 = readMatrix(sc, "B");
                        System.out.println("Результат деления:");
                        m1.divide(m2).print();
                        break;
                    }
                    case 4: {
                        Matrix m = readMatrix(sc, "A");
                        System.out.println("Транспонированная матрица:");
                        m.transpose().print();
                        break;
                    }
                    case 5: {
                        Matrix m = readMatrix(sc, "A");
                        System.out.println("Определитель матрицы: " + m.determinant());
                        break;
                    }
                    default:
                        System.out.println("Неверный пункт меню");
                }
            } catch (Exception e) {
                System.out.println("\nОшибка: " + e.getMessage());
                sc.nextLine();
            }
        }
    }
}