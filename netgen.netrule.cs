namespace netgen
{

	public class netrule
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public netrule()
		{
		  name = new char[1];
		  name[0] = (char)0;
		  quality = 0;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public void Dispose()
		{
		  name = null;
		  for (int i = 0; i < oldutofreearea_i.Size(); i++)
		  {
			oldutofreearea_i[i] = null;
		  }
		  for (int i = 0; i < freezone_i.Size(); i++)
		  {
			freezone_i[i] = null;
		  }
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public void SetFreeZoneTransformation(Vector devp, int tolclass)
		{
		  double lam1 = 1.0 / tolclass;
		  double lam2 = 1.0 - lam1;
        
		  double[] mem1 = new double[100];
		  double[] mem2 = new double[100];
		  double[] mem3 = new double[100];
        
		  int vs = oldutofreearea.Height();
		  FlatVector devfree = new FlatVector(vs, ref mem1);
        
		  int fzs = freezone.Size();
		  transfreezone.SetSize(fzs);
        
		  if (tolclass <= oldutofreearea_i.Size())
		  {
			  oldutofreearea_i[tolclass - 1].Mult(devp, devfree.functorMethod);
        
			  Array<Point2d> fzi = *freezone_i[tolclass - 1];
			  for (int i = 0; i < fzs; i++)
			  {
			  transfreezone[i].X() = fzi[i].X() + devfree.functorMethod[2 * i];
			  transfreezone[i].Y() = fzi[i].Y() + devfree.functorMethod[2 * i + 1];
			  }
		  }
		  else
		  {
			  FlatVector devfree1 = new FlatVector(vs, ref mem2);
			  FlatVector devfree2 = new FlatVector(vs, ref mem3);
        
			  oldutofreearea.Mult(devp, devfree1.functorMethod);
			  oldutofreearealimit.Mult(devp, devfree2.functorMethod);
			  devfree.Set2.functorMethod(lam1, devfree1.functorMethod, lam2, devfree2.functorMethod);
        
			  for (int i = 0; i < fzs; i++)
			  {
			  transfreezone[i].X() = lam1 * freezone[i].X() + lam2 * freezonelimit[i].X() + devfree.functorMethod[2 * i];
			  transfreezone[i].Y() = lam1 * freezone[i].Y() + lam2 * freezonelimit[i].Y() + devfree.functorMethod[2 * i + 1];
			  }
		  }
        
        
		  if (fzs > 0)
		  {
			  fzmaxx = fzminx = transfreezone[0].X();
			  fzmaxy = fzminy = transfreezone[0].Y();
		  }
        
		  for (int i = 1; i < fzs; i++)
		  {
			  if (transfreezone[i].X() > fzmaxx)
			  {
				  fzmaxx = transfreezone[i].X();
			  }
			  if (transfreezone[i].X() < fzminx)
			  {
				  fzminx = transfreezone[i].X();
			  }
			  if (transfreezone[i].Y() > fzmaxy)
			  {
				  fzmaxy = transfreezone[i].Y();
			  }
			  if (transfreezone[i].Y() < fzminy)
			  {
				  fzminy = transfreezone[i].Y();
			  }
		  }
        
		  for (int i = 0; i < fzs; i++)
		  {
			  Point2d p1 = transfreezone[i];
			  Point2d p2 = transfreezone[(i + 1) % fzs];
        
			  Vec2d vn = new Vec2d(p2.Y() - p1.Y(), p1.X() - p2.X());
        
			  double len2 = vn.Length2();
        
			  if (len2 < 1e-10)
			  {
			  freesetinequ(i, 0) = 0;
			  freesetinequ(i, 1) = 0;
			  freesetinequ(i, 2) = -1;
			  }
			  else
			  {
			  vn /= ngsimd.GlobalMembers.sqrt(len2); // scaling necessary ?
        
			  freesetinequ(i,0) = vn.X();
			  freesetinequ(i,1) = vn.Y();
			  freesetinequ(i,2) = -(p1.X() * vn.X() + p1.Y() * vn.Y());
			  }
		  }
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public int IsLineInFreeZone2(Point2d p1, Point2d p2)
		{
		  if ((p1.X() > fzmaxx && p2.X() > fzmaxx) || (p1.X() < fzminx && p2.X() < fzminx) || (p1.Y() > fzmaxy && p2.Y() > fzmaxy) || (p1.Y() < fzminy && p2.Y() < fzminy))
		  {
			  return 0;
		  }
        
		  for (int i = 1; i <= transfreezone.Size(); i++)
		  {
			  if (freesetinequ.Get(i, 1) * p1.X() + freesetinequ.Get(i, 2) * p1.Y() + freesetinequ.Get(i, 3) > -1e-8 && freesetinequ.Get(i, 1) * p2.X() + freesetinequ.Get(i, 2) * p2.Y() + freesetinequ.Get(i, 3) > -1e-8)
			  {
				  return 0;
			  }
		  }
        
		  double nx = (p2.Y() - p1.Y());
		  double ny = -(p2.X() - p1.X());
		  double nl = ngsimd.GlobalMembers.sqrt(nx * nx + ny * ny);
		  if (nl > 1e-8)
		  {
			  nx /= nl;
			  ny /= nl;
			  double c = - (p1.X() * nx + p1.Y() * ny);
        
			  bool allleft = true;
			  bool allright = true;
        
			  for (int i = 1; i <= transfreezone.Size(); i++)
			  {
			  bool left = transfreezone.Get(i).X() * nx + transfreezone.Get(i).Y() * ny + c < 1e-7;
				  bool right = transfreezone.Get(i).X() * nx + transfreezone.Get(i).Y() * ny + c > -1e-7;
			  if (!left)
			  {
				  allleft = false;
			  }
			  if (!right)
			  {
				  allright = false;
			  }
			  }
			  if (allleft || allright)
			  {
				  return false;
			  }
		  }
        
		  return true;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public int ConvexFreeZone()
		{
		  int n = transfreezone.Size();
		  for (int i = 1; i <= n; i++)
		  {
			  bool counterclockwise = CCW(transfreezone.Get(i), transfreezone.Get(i % n + 1), transfreezone.Get((i + 1) % n + 1), 1e-7);
			  //(*testout) << "ccw " << counterclockwise << endl << " p1 " << transfreezone.Get(i) << " p2 " << transfreezone.Get(i % n + 1)
			  //		 << " p3 " << transfreezone.Get( (i+1) % n + 1 ) << endl;
			  if (!counterclockwise)
			  {
			return 0;
			  }
		  }
		  return 1;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public float CalcLineError(int li, Vec2d v)
		{
		  float dx = v.X() - linevecs.Get(li).X();
		  float dy = v.Y() - linevecs.Get(li).Y();
        
		  threefloat ltf = linetolerances.Get(li);
		  return ltf.f1 * dx * dx + ltf.f2 * dx * dy + ltf.f3 * dy * dy;
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public void LoadRule(istream ist)
		{
		  string buf = new string(new char[256]);
		  char ch;
		  Point2d p = new Point2d();
		  INDEX_2 lin = new INDEX_2();
		  int i;
		  int j;
		  DenseMatrix tempoldutonewu = new DenseMatrix(20, 20);
		  DenseMatrix tempoldutofreearea = new DenseMatrix(20, 20);
		  DenseMatrix tempoldutofreearealimit = new DenseMatrix(20, 20);
        
		  tempoldutonewu = 0;
		  tempoldutofreearea = 0;
		  tempoldutofreearealimit = 0;
        
		  noldp = 0;
		  noldl = 0;
        
		  ist.get(buf, sizeof(char), '"');
		  ist.get(ch);
		  ist.get(buf, sizeof(char), '"');
		  ist.get(ch);
        
		  // if(name != NULL) 
		  name = null;
		  name = new char[buf.Length + 1];
		  name = buf;
		  //(*testout) << "name " << name << endl;
		  //  (*mycout) << "Rule " << name << " found." << endl;
        
		  do
		  {
			  ist >> buf;
        
			  //(*testout) << "buf " << buf << endl;
        
			  if (string.Compare(buf, "quality") == 0)
        
			  {
			  ist >> quality;
			  }
        
			  else if (string.Compare(buf, "mappoints") == 0)
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  ist >> p.X();
				  ist >> ch; // ','
				  ist >> p.Y();
				  ist >> ch; // ')'
        
				  points.Append(p);
				  noldp++;
        
				  tolerances.SetSize(noldp);
				  tolerances.Elem(noldp).f1 = 1.0;
				  tolerances.Elem(noldp).f2 = 0;
				  tolerances.Elem(noldp).f3 = 1.0;
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  if (ch == '{')
				  {
					  ist >> tolerances.Elem(noldp).f1;
					  ist >> ch; // ','
					  ist >> tolerances.Elem(noldp).f2;
					  ist >> ch; // ','
					  ist >> tolerances.Elem(noldp).f3;
					  ist >> ch; // '}'
				  }
				  else if (ch == 'd')
				  {
					  //            delpoints.Append (noldp);
					  ist >> ch; // 'e'
					  ist >> ch; // 'l'
				  }
        
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  ist.putback(ch);
			  }
        
        
			  else if (string.Compare(buf, "maplines") == 0)
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  ist >> lin.I1();
				  ist >> ch; // ','
				  ist >> lin.I2();
				  ist >> ch; // ')'
        
        
				  //(*testout) << "read line " << lin.I1() << " " << lin.I2() << endl;
				  lines.Append(lin);
				  linevecs.Append(points.Get(lin.I2()) - points.Get(lin.I1()));
				  noldl++;
				  linetolerances.SetSize(noldl);
				  linetolerances.Elem(noldl).f1 = 0;
				  linetolerances.Elem(noldl).f2 = 0;
				  linetolerances.Elem(noldl).f3 = 0;
        
				  //(*testout) << "mapl1" << endl; 
				  ist >> ch;
				  while (ch != ';')
				  {
				  //(*testout) << "working on character \""<<ch<<"\""<< endl;
				  if (ch == '{')
				  {
					  ist >> linetolerances.Elem(noldl).f1;
					  ist >> ch; // ','
					  ist >> linetolerances.Elem(noldl).f2;
					  ist >> ch; // ','
					  ist >> linetolerances.Elem(noldl).f3;
					  ist >> ch; // '}'
				  }
				  else if (ch == 'd')
				  {
					  dellines.Append(noldl);
					  ist >> ch; // 'e'
					  ist >> ch; // 'l'
					  //(*testout) << "read del" << endl;
				  }
        
				  ist >> ch;
				  //(*testout) << "read character \""<<ch<<"\""<< endl;
				  }
        
				  ist >> ch;
				  //(*testout) << "read next character \""<<ch<<"\""<< endl;
			  }
        
        
			  ist.putback(ch);
			  }
        
			  else if (string.Compare(buf, "newpoints") == 0)
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  ist >> p.X();
				  ist >> ch; // ','
				  ist >> p.Y();
				  ist >> ch; // ')'
        
				  points.Append(p);
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  if (ch == '{')
				  {
					  LoadMatrixLine(ist, tempoldutonewu.functorMethod, 2 * (points.Size() - noldp) - 1);
        
					  ist >> ch; // '{'
					  LoadMatrixLine(ist, tempoldutonewu.functorMethod, 2 * (points.Size() - noldp));
				  }
        
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  ist.putback(ch);
			  }
        
			  else if (string.Compare(buf, "newlines") == 0)
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  ist >> lin.I1();
				  ist >> ch; // ','
				  ist >> lin.I2();
				  ist >> ch; // ')'
        
				  lines.Append(lin);
				  linevecs.Append(points.Get(lin.I2()) - points.Get(lin.I1()));
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  ist.putback(ch);
			  }
        
			  else if (string.Compare(buf, "freearea") == 0)
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  ist >> p.X();
				  ist >> ch; // ','
				  ist >> p.Y();
				  ist >> ch; // ')'
        
				  freezone.Append(p);
				  freezonelimit.Append(p);
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  if (ch == '{')
				  {
					  LoadMatrixLine(ist, tempoldutofreearea.functorMethod, 2 * freezone.Size() - 1);
        
					  ist >> ch; // '{'
					  LoadMatrixLine(ist, tempoldutofreearea.functorMethod, 2 * freezone.Size());
				  }
        
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  for (i = 1; i <= tempoldutofreearealimit.Height(); i++)
			  {
				for (j = 1; j <= tempoldutofreearealimit.Width(); j++)
				{
				  tempoldutofreearealimit.Elem(i, j) = tempoldutofreearea.Elem(i, j);
				}
			  }
        
        
			  ist.putback(ch);
			  }
			  else if (string.Compare(buf, "freearea2") == 0)
			  {
			  ist >> ch;
			  int freepi = 0;
			  tempoldutofreearealimit = 0;
        
			  while (ch == '(')
			  {
				  freepi++;
        
				  ist >> p.X();
				  ist >> ch; // ','
				  ist >> p.Y();
				  ist >> ch; // ')'
        
				  freezonelimit.Elem(freepi) = p;
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  if (ch == '{')
				  {
					  LoadMatrixLine(ist, tempoldutofreearealimit.functorMethod, 2 * freepi - 1);
        
					  ist >> ch; // '{'
					  LoadMatrixLine(ist, tempoldutofreearealimit.functorMethod, 2 * freepi);
				  }
        
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  ist.putback(ch);
			  }
        
			  else if (string.Compare(buf, "elements") == 0)
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  elements.Append(new Element2d(ELEMENT_TYPE.TRIG));
        
				  ist >> elements.Last().PNum(1);
				  ist >> ch; // ','
        
				  if (ch == DefineConstants.COMMASIGN)
				  {
				  ist >> elements.Last().PNum(2);
				  ist >> ch; // ','
				  }
				  if (ch == DefineConstants.COMMASIGN)
				  {
				  ist >> elements.Last().PNum(3);
				  ist >> ch; // ','
				  }
				  if (ch == DefineConstants.COMMASIGN)
				  {
				  elements.Last().SetType(ELEMENT_TYPE.QUAD);
				  ist >> elements.Last().PNum(4);
				  ist >> ch; // ','
        
				  // const Element2d & el = elements.Last();
				  /*
				  orientations.Append (threeint(el.PNum(1), el.PNum(2), el.PNum(3)));
				  orientations.Append (threeint(el.PNum(2), el.PNum(3), el.PNum(4)));
				  orientations.Append (threeint(el.PNum(3), el.PNum(4), el.PNum(1)));
				  orientations.Append (threeint(el.PNum(4), el.PNum(1), el.PNum(2)));
				  */
				  }
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  ist.putback(ch);
			  }
        
			  else if (string.Compare(buf, "orientations") == 0)
        
			  {
			  ist >> ch;
        
			  while (ch == '(')
			  {
				  //        threeint a = threeint();
				  orientations.Append(new threeint());
        
				  ist >> orientations.Last().i1;
				  ist >> ch; // ','
				  ist >> orientations.Last().i2;
				  ist >> ch; // ','
				  ist >> orientations.Last().i3;
				  ist >> ch; // ','
        
				  ist >> ch;
				  while (ch != ';')
				  {
				  ist >> ch;
				  }
        
				  ist >> ch;
			  }
        
			  ist.putback(ch);
			  }
        
			  else if (string.Compare(buf, "endrule") != 0)
			  {
			  PrintSysError("Parser error, unknown token ", buf);
			  }
		  } while (!ist.eof() && string.Compare(buf, "endrule") != 0);
        
		  oldutonewu.SetSize(2 * (points.Size() - noldp), 2 * noldp);
		  oldutofreearea.SetSize(2 * freezone.Size(), 2 * noldp);
		  oldutofreearealimit.SetSize(2 * freezone.Size(), 2 * noldp);
        
		  for (i = 1; i <= oldutonewu.Height(); i++)
		  {
			for (j = 1; j <= oldutonewu.Width(); j++)
			{
			  oldutonewu.Elem(i, j) = tempoldutonewu.Elem(i, j);
			}
		  }
        
		  for (i = 1; i <= oldutofreearea.Height(); i++)
		  {
			for (j = 1; j <= oldutofreearea.Width(); j++)
			{
			  oldutofreearea.Elem(i, j) = tempoldutofreearea.Elem(i, j);
			}
		  }
        
		  for (i = 1; i <= oldutofreearea.Height(); i++)
		  {
			for (j = 1; j <= oldutofreearea.Width(); j++)
			{
			  oldutofreearealimit.Elem(i, j) = tempoldutofreearealimit.Elem(i, j);
			}
		  }
        
		  freesetinequ.SetSize(freezone.Size());
        
        
		  {
			char ok;
			int minn;
			Array<int> pnearness = new Array<int>(noldp);
        
			for (i = 1; i <= pnearness.Size(); i++)
			{
			  pnearness.Elem(i) = 1000;
			}
        
			for (j = 1; j <= 2; j++)
			{
			  pnearness.Elem(GetPointNr(1, j)) = 0;
			}
        
			do
			{
			ok = 1;
        
			for (i = 1; i <= noldl; i++)
			{
				minn = 1000;
				for (j = 1; j <= 2; j++)
				{
				  minn = min2(minn, pnearness.Get(GetPointNr(i, j)));
				}
        
				for (j = 1; j <= 2; j++)
				{
				  if (pnearness.Get(GetPointNr(i, j)) > minn + 1)
				  {
				  ok = 0;
				  pnearness.Elem(GetPointNr(i, j)) = minn + 1;
				  }
				}
			}
			} while (!ok);
        
			lnearness.SetSize(noldl);
        
			for (i = 1; i <= noldl; i++)
			{
			lnearness.Elem(i) = 0;
			for (j = 1; j <= 2; j++)
			{
			  lnearness.Elem(i) += pnearness.Get(GetPointNr(i, j));
			}
			}
		  }
        
		  oldutofreearea_i.SetSize(10);
		  freezone_i.SetSize(10);
        
		  for (i = 0; i < oldutofreearea_i.Size(); i++)
		  {
			  double lam1 = 1.0 / (i + 1);
        
			  oldutofreearea_i[i] = new DenseMatrix(oldutofreearea.Height(), oldutofreearea.Width());
			  DenseMatrix mati = *oldutofreearea_i[i];
			  for (j = 0; j < oldutofreearea.Height(); j++)
			  {
			for (int k = 0; k < oldutofreearea.Width(); k++)
			{
			  mati.functorMethod(j,k) = lam1 * oldutofreearea(j,k) + (1 - lam1) * oldutofreearealimit(j,k);
			}
			  }
        
			  freezone_i[i] = new Array<Point2d> (freezone.Size());
			  Array<Point2d> fzi = *freezone_i[i];
			  for (int j = 0; j < freezone.Size(); j++)
			  {
			fzi[j] = freezonelimit[j] + lam1 * (freezone[j] - freezonelimit[j]);
			  }
		  }
		}
	}
}