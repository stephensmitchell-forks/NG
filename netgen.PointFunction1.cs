namespace netgen
{

	public class PointFunction1
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public PointFunction1(Mesh.T_POINTS apoints, Array<INDEX_3> afaces, MeshingParameters amp, double ah)
		  {
			  this.points = apoints;
			  this.faces = afaces;
			  this.mp = amp;
			h = ah;
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public double Func(Vector vp)
		  {
			double badness = 0;
			Point < 3> pp(vp(0), vp(1), vp(2));
        
			for (int j = 0; j < faces.Size(); j++)
			{
			INDEX_3 el = faces[j];
        
			double bad = CalcTetBadness(points[new PointIndex(el.I1())], points[new PointIndex(el.I3())], points[new PointIndex(el.I2())], pp, 0, mp);
			badness += bad;
			}
        
			return badness;
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public double FuncDeriv(Vector x, Vector dir, ref double deriv)
		  {
			VectorMem < 3> hx;
			const double eps = 1e-6;
        
			double dirlen = dir.L2Norm();
			if (dirlen < 1e-14)
			{
			deriv = 0;
			return Func(x);
			}
        
			hx.Set(1, x);
			hx.Add(eps * h / dirlen, dir);
			double fr = Func(hx);
			hx.Set(1, x);
			hx.Add(-eps * h / dirlen, dir);
			double fl = Func(hx);
        
			deriv = (fr - fl) / (2 * eps * h) * dirlen;
        
			return Func(x);
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public double FuncGrad(Vector x, Vector g)
		  {
			VectorMem < 3> hx;
			double eps = 1e-6;
        
			hx = x;
			for (int i = 0; i < 3; i++)
			{
			hx(i) = x(i) + eps * h;
			double fr = Func(hx);
			hx(i) = x(i) - eps * h;
			double fl = Func(hx);
			hx(i) = x(i);
        
			g(i) = (fr - fl) / (2 * eps * h);
			}
        
			return Func(x);
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public double GradStopping(Vector x)
		  {
			double f = Func(x);
			return 1e-8 * f * f;
		  }
	}
}