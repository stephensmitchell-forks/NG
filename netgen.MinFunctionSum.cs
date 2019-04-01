namespace netgen
{

	public class MinFunctionSum
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public double Func(Vector x)
		  {
			double retval = 0;
			for (int i = 0; i < functions.Size(); i++)
			{
			  retval += functions[i].Func(x);
			}
        
			return retval;
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void Grad(Vector x, ref Vector g)
		  {
			g = 0.0;
			VectorMem < 3> gi;
			for (int i = 0; i < functions.Size(); i++)
			{
			functions[i].Grad(x,gi);
			for (int j = 0; j < g.Size(); j++)
			{
			  g[j] += gi[j];
			}
			}
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public double FuncGrad(Vector x, ref Vector g)
		  {
			double retval = 0;
			g = 0.0;
			VectorMem < 3> gi;
			for (int i = 0; i < functions.Size(); i++)
			{
			retval += functions[i].FuncGrad(x,gi);
			for (int j = 0; j < g.Size(); j++)
			{
			  g[j] += gi[j];
			}
			}
			return retval;
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public double FuncDeriv(Vector x, Vector dir, ref double deriv)
		  {
			double retval = 0;
			deriv = 0.0;
			double derivi;
			for (int i = 0; i < functions.Size(); i++)
			{
			retval += functions[i].FuncDeriv(x,dir,derivi);
			deriv += derivi;
			}
			return retval;
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public double GradStopping(Vector x)
		  {
			double minfs = 0;
			double mini;
			for (int i = 0; i < functions.Size(); i++)
			{
			mini = functions[i].GradStopping(x);
			if (i == 0 || mini < minfs)
			{
			  minfs = mini;
			}
			}
			return minfs;
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void AddFunction(MinFunction fun)
		  {
			functions.Append(fun);
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public MinFunction Function(int i)
		  {
			return functions[i];
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public MinFunction Function(int i)
		  {
			return functions[i];
		  }
	}
}