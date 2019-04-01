namespace netgen
{

	public class Element2d
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GetTransformation(int ip, Array<Point2d> points, DenseMatrix trans)
		  {
			int np = GetNP();
			DenseMatrix pmat = new DenseMatrix(2, np);
			DenseMatrix dshape = new DenseMatrix(2, np);
			pmat.SetSize(2, np);
			dshape.SetSize(2, np);
        
			Point2d p = new Point2d();
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
		  public void GetShape(Point2d p, Vector shape)
		  {
			if (shape.Size() != GetNP())
			{
				cerr << "Element::GetShape: Length not fitting" << "\n";
				return;
			}
        
			switch (typ)
			{
			  case ELEMENT_TYPE.TRIG:
				shape(0) = 1 - p.X() - p.Y();
				shape(1) = p.X();
				shape(2) = p.Y();
				break;
			  case ELEMENT_TYPE.QUAD:
				shape(0) = (1 - p.X()) * (1 - p.Y());
				shape(1) = p.X() * (1 - p.Y());
				shape(2) = p.X() * p.Y();
				shape(3) = (1 - p.X()) * p.Y();
				break;
			  default:
				PrintSysError("Element2d::GetShape, illegal type ", (int)typ);
				break;
			}
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GetShapeNew(Point < 2> p, FlatVector shape)
		  {
			switch (typ)
			{
			  case ELEMENT_TYPE.TRIG:
			  {
				  shape.functorMethod(0) = p.functorMethod(0);
				  shape.functorMethod(1) = p.functorMethod(1);
				  shape.functorMethod(2) = 1 - p.functorMethod(0) - p.functorMethod(1);
				  break;
			  }
        
			  case ELEMENT_TYPE.QUAD:
			  {
				  shape.functorMethod(0) = (1 - p.functorMethod(0)) * (1 - p.functorMethod(1));
				  shape.functorMethod(1) = p.functorMethod(0) * (1 - p.functorMethod(1));
				  shape.functorMethod(2) = p.functorMethod(0) * p.functorMethod(1);
				  shape.functorMethod(3) = (1 - p.functorMethod(0)) * p.functorMethod(1);
				  break;
			  }
        
			  default:
				throw new Exception("illegal element type in GetShapeNew");
			}
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GetDShape(Point2d p, DenseMatrix dshape)
		  {
		#if DEBUG
			if (dshape.Height() != 2 || dshape.Width() != np)
			{
				PrintSysError("Element::DShape: Sizes don't fit");
				return;
			}
		#endif
        
			switch (typ)
			{
			  case ELEMENT_TYPE.TRIG:
				dshape.Elem(1, 1) = -1;
				dshape.Elem(1, 2) = 1;
				dshape.Elem(1, 3) = 0;
				dshape.Elem(2, 1) = -1;
				dshape.Elem(2, 2) = 0;
				dshape.Elem(2, 3) = 1;
				break;
			  case ELEMENT_TYPE.QUAD:
				dshape.Elem(1, 1) = -(1 - p.Y());
				dshape.Elem(1, 2) = (1 - p.Y());
				dshape.Elem(1, 3) = p.Y();
				dshape.Elem(1, 4) = -p.Y();
				dshape.Elem(2, 1) = -(1 - p.X());
				dshape.Elem(2, 2) = -p.X();
				dshape.Elem(2, 3) = p.X();
				dshape.Elem(2, 4) = (1 - p.X());
				break;
        
			  default:
				PrintSysError("Element2d::GetDShape, illegal type ", (int)typ);
				break;
			}
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GetDShapeNew<T>(Point<2,T> p, ref MatrixFixWidth<2,T> dshape)
		  {
			switch (typ)
			{
			  case ELEMENT_TYPE.TRIG:
			  {
				  dshape = T(0.0);
				  dshape(0,0) = 1;
				  dshape(1,1) = 1;
				  dshape(2,0) = -1;
				  dshape(2,1) = -1;
				  break;
			  }
			  case ELEMENT_TYPE.QUAD:
			  {
				  dshape(0,0) = -(1 - p(1));
				  dshape(0,1) = -(1 - p(0));
        
				  dshape(1,0) = (1 - p(1));
				  dshape(1,1) = -p(0);
        
				  dshape(2,0) = p(1);
				  dshape(2,1) = p(0);
        
				  dshape(3,0) = -p(1);
				  dshape(3,1) = (1 - p(0));
				  break;
			  }
			  default:
				throw new Exception("illegal element type in GetDShapeNew");
			}
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GetPointMatrix(Array<Point2d> points, DenseMatrix pmat)
		  {
			int np = GetNP();
        
		#if DEBUG
			if (pmat.Width() != np || pmat.Height() != 2)
			{
				cerr << "Element::GetPointMatrix: sizes don't fit" << "\n";
				return;
			}
		#endif
        
			for (int i = 1; i <= np; i++)
			{
				Point2d p = points.Get(PNum(i));
				pmat.Elem(1, i) = p.X();
				pmat.Elem(2, i) = p.Y();
			}
		  }
	}
}