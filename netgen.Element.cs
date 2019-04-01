namespace netgen
{

	public class Element
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GetNodesLocal(Array<Point3d> points)
		  {
			double[][] tetpoints =
			{
				new double[] {0, 0, 0},
				new double[] {1, 0, 0},
				new double[] {0, 1, 0},
				new double[] {0, 0, 1}
			};
        
			double[][] prismpoints =
			{
				new double[] {0, 0, 0},
				new double[] {1, 0, 0},
				new double[] {0, 1, 0},
				new double[] {0, 0, 1},
				new double[] {1, 0, 1},
				new double[] {0, 1, 1}
			};
        
			double[][] pyramidpoints =
			{
				new double[] {0, 0, 0},
				new double[] {1, 0, 0},
				new double[] {1, 1, 0},
				new double[] {0, 1, 0},
				new double[] {0, 0, 1}
			};
        
			double[][] tet10points =
			{
				new double[] {0, 0, 0},
				new double[] {1, 0, 0},
				new double[] {0, 1, 0},
				new double[] {0, 0, 1},
				new double[] {0.5, 0, 0},
				new double[] {0, 0.5, 0},
				new double[] {0, 0, 0.5},
				new double[] {0.5, 0.5, 0},
				new double[] {0.5, 0, 0.5},
				new double[] {0, 0.5, 0.5}
			};
        
			double[][] hexpoints =
			{
				new double[] {0, 0, 0},
				new double[] {1, 0, 0},
				new double[] {1, 1, 0},
				new double[] {0, 1, 0},
				new double[] {0, 0, 1},
				new double[] {1, 0, 1},
				new double[] {1, 1, 1},
				new double[] {0, 1, 1}
			};
        
			int np;
			int i;
			double[] pp = new double[3];
			switch (GetType())
			{
			  case ELEMENT_TYPE.TET:
			  {
				  np = 4;
				  pp = tetpoints;
				  break;
			  }
			  case ELEMENT_TYPE.PRISM:
			  case ELEMENT_TYPE.PRISM12:
			  {
				  np = 6;
				  pp = prismpoints;
				  break;
			  }
			  case ELEMENT_TYPE.TET10:
			  {
				  np = 10;
				  pp = tet10points;
				  break;
			  }
			  case ELEMENT_TYPE.PYRAMID:
			  {
				  np = 5;
				  pp = pyramidpoints;
				  break;
			  }
			  case ELEMENT_TYPE.HEX:
			  {
				  np = 8;
				  pp = hexpoints;
				  break;
			  }
			  default:
			  {
				  Console.Write("GetNodesLocal not impelemented for element ");
				  Console.Write(GetType());
				  Console.Write("\n");
				  np = 0;
			  }
				break;
			}
        
			points.SetSize(0);
			for (i = 0; i < np; i++)
			{
			  points.Append(new Point3d(pp[i][0], pp[i][1], pp[i][2]));
			}
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GetTransformation(int ip, T_POINTS points, DenseMatrix trans)
		  {
			int np = GetNP();
			DenseMatrix pmat = new DenseMatrix(3, np);
			DenseMatrix dshape = new DenseMatrix(3, np);
			pmat.SetSize(3, np);
			dshape.SetSize(3, np);
        
			Point < 3> p;
			double w;
        
			GetPointMatrix(points, pmat.functorMethod);
			GetIntegrationPoint(ip, p, w);
			GetDShape(p, dshape.functorMethod);
        
			CalcABt(pmat.functorMethod, dshape.functorMethod, trans.functorMethod);
        
			/*
			  (*testout) << "p = " << p  << endl
			  << "pmat = " << pmat << endl
			  << "dshape = " << dshape << endl
			  << "tans = " << trans << endl;
			*/
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GetShape(Point < 3> hp, Vector shape)
		  {
			Point3d p = hp.functorMethod;
        
			if (shape.Size() != GetNP())
			{
				cerr << "Element::GetShape: Length not fitting" << "\n";
				return;
			}
        
			switch (typ)
			{
			  case ELEMENT_TYPE.TET:
			  {
				  shape(0) = 1 - p.X() - p.Y() - p.Z();
				  shape(1) = p.X();
				  shape(2) = p.Y();
				  shape(3) = p.Z();
				  break;
			  }
			  case ELEMENT_TYPE.TET10:
			  {
				  double lam1 = 1 - p.X() - p.Y() - p.Z();
				  double lam2 = p.X();
				  double lam3 = p.Y();
				  double lam4 = p.Z();
        
				  shape(4) = 4 * lam1 * lam2;
				  shape(5) = 4 * lam1 * lam3;
				  shape(6) = 4 * lam1 * lam4;
				  shape(7) = 4 * lam2 * lam3;
				  shape(8) = 4 * lam2 * lam4;
				  shape(9) = 4 * lam3 * lam4;
        
				  shape(0) = lam1 - 0.5 * (shape(4) + shape(5) + shape(6));
				  shape(1) = lam2 - 0.5 * (shape(4) + shape(7) + shape(8));
				  shape(2) = lam3 - 0.5 * (shape(5) + shape(7) + shape(9));
				  shape(3) = lam4 - 0.5 * (shape(6) + shape(8) + shape(9));
				  break;
			  }
        
			  case ELEMENT_TYPE.PRISM:
			  {
				  Point < 3> hp = p;
				  shape(0) = hp.functorMethod(0) * (1 - hp.functorMethod(2));
				  shape(1) = hp.functorMethod(1) * (1 - hp.functorMethod(2));
				  shape(2) = (1 - hp.functorMethod(0) - hp.functorMethod(1)) * (1 - hp.functorMethod(2));
				  shape(3) = hp.functorMethod(0) * hp.functorMethod(2);
				  shape(4) = hp.functorMethod(1) * hp.functorMethod(2);
				  shape(5) = (1 - hp.functorMethod(0) - hp.functorMethod(1)) * hp.functorMethod(2);
				  break;
			  }
			  case ELEMENT_TYPE.HEX:
			  {
				  Point < 3> hp = p;
				  shape(0) = (1 - hp.functorMethod(0)) * (1 - hp.functorMethod(1)) * (1 - hp.functorMethod(2));
				  shape(1) = .functorMethod(hp.functorMethod(0)) * (1 - hp.functorMethod(1)) * (1 - hp.functorMethod(2));
				  shape(2) = .functorMethod(hp.functorMethod(0)) * .functorMethod(hp.functorMethod(1)) * (1 - hp.functorMethod(2));
				  shape(3) = (1 - hp.functorMethod(0)) * .functorMethod(hp.functorMethod(1)) * (1 - hp.functorMethod(2));
				  shape(4) = (1 - hp.functorMethod(0)) * (1 - hp.functorMethod(1)) * .functorMethod(hp.functorMethod(2));
				  shape(5) = .functorMethod(hp.functorMethod(0)) * (1 - hp.functorMethod(1)) * .functorMethod(hp.functorMethod(2));
				  shape(6) = .functorMethod(hp.functorMethod(0)) * .functorMethod(hp.functorMethod(1)) * .functorMethod(hp.functorMethod(2));
				  shape(7) = (1 - hp.functorMethod(0)) * .functorMethod(hp.functorMethod(1)) * .functorMethod(hp.functorMethod(2));
				  break;
			  }
			  default:
				throw new Exception("Element :: GetShape not implemented for that element");
			}
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GetDShape(Point < 3> hp, DenseMatrix dshape)
		  {
			Point3d p = hp.functorMethod;
        
			int np = GetNP();
			if (dshape.Height() != 3 || dshape.Width() != np)
			{
				cerr << "Element::DShape: Sizes don't fit" << "\n";
				return;
			}
        
			double eps = 1e-6;
			Vector shaper = new Vector(np);
			Vector shapel = new Vector(np);
        
			for (int i = 1; i <= 3; i++)
			{
				Point3d pr = new Point3d(p);
				Point3d pl = new Point3d(p);
				pr.X(i) += eps;
				pl.X(i) -= eps;
        
				GetShape(pr, shaper);
				GetShape(pl, shapel);
				for (int j = 0; j < np; j++)
				{
				  dshape.functorMethod(i - 1, j) = (shaper(j) - shapel(j)) / (2 * eps);
				}
			}
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GetDShapeNew<T>(Point<3,T> p, ref MatrixFixWidth<3,T> dshape)
		  {
			switch (typ)
			{
			  case ELEMENT_TYPE.TET:
			  {
				  dshape = T(0.0);
				  dshape(0,0) = 1;
				  dshape(1,1) = 1;
				  dshape(2,2) = 1;
				  dshape(3,0) = -1;
				  dshape(3,1) = -1;
				  dshape(3,2) = -1;
				  break;
			  }
			  case ELEMENT_TYPE.PRISM:
			  {
				  dshape = T(0.0);
				  dshape(0,0) = 1 - p(2);
				  dshape(0,2) = -p(0);
				  dshape(1,1) = 1 - p(2);
				  dshape(1,2) = -p(1);
				  dshape(2,0) = -(1 - p(2));
				  dshape(2,1) = -(1 - p(2));
				  dshape(2,2) = -(1 - p(0) - p(1));
        
				  dshape(3,0) = p(2);
				  dshape(3,2) = p(0);
				  dshape(4,1) = p(2);
				  dshape(4,2) = p(1);
				  dshape(5,0) = -p(2);
				  dshape(5,1) = -p(2);
				  dshape(5,2) = 1 - p(0) - p(1);
				  break;
			  }
        
			  default:
			  {
				  int np = GetNP();
				  double eps = 1e-6;
				  ArrayMem<T,100> mem = new ArrayMem<T,100>((uint)(2 * np));
				  TFlatVector<T> shaper = new TFlatVector<T>(np, mem[0]);
				  TFlatVector<T> shapel = new TFlatVector<T>(np, mem[np]);
				  // Vector shaper(np), shapel(np);
        
				  for (int i = 0; i < 3; i++)
				  {
					  Point<3,T> pr = new Point<3,T>(p);
					  Point<3,T> pl = new Point<3,T>(p);
					  pr(i) += eps;
					  pl(i) -= eps;
        
					  GetShapeNew(pr, shaper.functorMethod);
					  GetShapeNew(pl, shapel.functorMethod);
					  for (int j = 0; j < np; j++)
					  {
						dshape(j, i) = (shaper.functorMethod(j) - shapel.functorMethod(j)) / (2 * eps);
					  }
				  }
			  }
				break;
			}
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GetPointMatrix(T_POINTS points, DenseMatrix pmat)
		  {
			int np = GetNP();
			for (int i = 1; i <= np; i++)
			{
				Point3d p = points.Get(PNum(i));
				pmat.Elem(1, i) = p.X();
				pmat.Elem(2, i) = p.Y();
				pmat.Elem(3, i) = p.Z();
			}
		  }
	}
}