namespace netgen
{

	public class JacobianPointFunction
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public JacobianPointFunction(Mesh.T_POINTS apoints, Array<Element, 0, uint> aelements)
		{
			this.points = apoints;
			this.elements = aelements;
			this.elementsonpoint = apoints.Size();
		  int i;
		  int j;
        
		  for (i = 1; i <= elements.Size(); i++)
		  {
			  for (j = 1; j <= elements.Get(i).NP(); j++)
			  {
			elementsonpoint.Add1(elements.Get(i).PNum(j), i);
			  }
		  }
        
		  onplane = false;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public void SetPointIndex(PointIndex aactpind)
		{
		  actpind = aactpind;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public double Func(Vector v)
		{
		  int j;
		  double badness = 0;
        
		  Point < 3> hp = points.Elem(actpind);
        
		  points.Elem(actpind) = hp + Vec < 3> (v(0), v(1), v(2));
        
		  if (onplane)
		  {
			points.Elem(actpind) -= (v(0) * nv(0) + v(1) * nv(1) + v(2) * nv(2)) * nv;
		  }
        
        
		  for (j = 1; j <= elementsonpoint.EntrySize(actpind); j++)
		  {
			  int eli = elementsonpoint.Get(actpind, j);
			  badness += elements.Get(eli).CalcJacobianBadness(points);
		  }
        
		  points.Elem(actpind) = hp;
        
		  return badness;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public double FuncGrad(Vector x, ref Vector g)
		{
		  int j;
		  int k;
		  int lpi;
		  double badness = 0; //, hbad;
        
		  Point < 3> hp = points.Elem(actpind);
		  points.Elem(actpind) = hp + Vec < 3> (x(0), x(1), x(2));
        
		  if (onplane)
		  {
			points.Elem(actpind) -= (x(0) * nv(0) + x(1) * nv(1) + x(2) * nv(2)) * nv;
		  }
        
		  Vec < 3> hderiv;
		  //Vec3d vdir;
		  g.SetSize(3);
		  g = 0;
        
		  for (j = 1; j <= elementsonpoint.EntrySize(actpind); j++)
		  {
			  int eli = elementsonpoint.Get(actpind, j);
			  Element el = elements.Get(eli);
        
			  lpi = 0;
			  for (k = 1; k <= el.GetNP(); k++)
			  {
			if (el.PNum(k) == actpind)
			{
			  lpi = k;
			}
			  }
			  if (lpi == 0)
			  {
				  cerr << "loc point not found" << "\n";
			  }
        
			  badness += elements.Get(eli).CalcJacobianBadnessGradient(points, lpi, hderiv);
        
			  for (k = 0; k < 3; k++)
			  {
			g(k) += hderiv(k);
			  }
        
			  /*
			  for (k = 1; k <= 3; k++)
			{
			  vdir = Vec3d(0,0,0);
			  vdir.X(k) = 1;
		
			  hbad = elements.Get(eli).
				CalcJacobianBadnessDirDeriv (points, lpi, vdir, hderiv);
			  //(*testout) << "hderiv " << k << ": " << hderiv << endl;
			  g.Elem(k) += hderiv;
			  if (k == 1)
				badness += hbad;
			}
			  */
		  }
        
		  if (onplane)
		  {
			  double scal = nv(0) * g(0) + nv(1) * g(1) + nv(2) * g(2);
			  g(0) -= scal * nv(0);
			  g(1) -= scal * nv(1);
			  g(2) -= scal * nv(2);
		  }
        
		  //(*testout) << "g = " << g << endl;
        
        
		  points.Elem(actpind) = hp;
        
		  return badness;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public double FuncDeriv(Vector x, Vector dir, ref double deriv)
		{
		  int j;
		  int k;
		  int lpi;
		  double badness = 0;
        
		  Point < 3> hp = points.Elem(actpind);
		  points.Elem(actpind) = Point < 3> (hp + new Vec3d(x(0), x(1), x(2)));
        
		  if (onplane)
		  {
			points.Elem(actpind) -= (new Vec3d(x(0), x(1), x(2)) * nv) * nv;
		  }
        
		  double hderiv;
		  deriv = 0;
		  Vec < 3> vdir(dir(0), dir(1), dir(2));
        
		  if (onplane)
		  {
			  double scal = vdir * nv;
			  vdir -= scal * nv;
		  }
        
		  for (j = 1; j <= elementsonpoint.EntrySize(actpind); j++)
		  {
			  int eli = elementsonpoint.Get(actpind, j);
			  Element el = elements.Get(eli);
        
			  lpi = 0;
			  for (k = 1; k <= el.GetNP(); k++)
			  {
			if (el.PNum(k) == actpind)
			{
			  lpi = k;
			}
			  }
			  if (lpi == 0)
			  {
				  cerr << "loc point not found" << "\n";
			  }
        
			  badness += elements.Get(eli).CalcJacobianBadnessDirDeriv(points, lpi, vdir, hderiv);
			  deriv += hderiv;
		  }
        
		  points.Elem(actpind) = hp;
        
		  return badness;
        
		}
	}
}