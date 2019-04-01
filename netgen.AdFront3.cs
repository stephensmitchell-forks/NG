namespace netgen
{

	public class AdFront3
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public void DeleteFace(int fi)
		{
		  nff--;
        
		  for (int i = 1; i <= faces.Get(fi).Face().GetNP(); i++)
		  {
			  PointIndex pi = faces.Get(fi).Face().PNum(i);
			  points[pi].RemoveFace();
			  if (!points[pi].Valid())
			  {
			delpointl.Append(pi);
			  }
		  }
        
		  MiniElement2d face = faces.Get(fi).Face();
		  Point3d p1 = points[face.PNum(1)].P();
		  Point3d p2 = points[face.PNum(2)].P();
		  Point3d p3 = points[face.PNum(3)].P();
        
		  vol -= 1.0 / 6.0 * (p1.X() + p2.X() + p3.X()) * ((p2.Y() - p1.Y()) * (p3.Z() - p1.Z()) - (p2.Z() - p1.Z()) * (p3.Y() - p1.Y()));
        
		  if (face.GetNP() == 4)
		  {
			  Point3d p4 = points[face.PNum(4)].P();
			  vol -= 1.0 / 6.0 * (p1.X() + p3.X() + p4.X()) * ((p3.Y() - p1.Y()) * (p4.Z() - p1.Z()) - (p3.Z() - p1.Z()) * (p4.Y() - p1.Y()));
        
			  nff4--;
		  }
        
		  faces.Elem(fi).Invalidate();
		}

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		public int GetLocals(int fstind, Array<Point3d, PointIndex.BASE> locpoints, Array<MiniElement2d> locfaces, Array<PointIndex, PointIndex.BASE> pindex, Array<int> findex, INDEX_2_HASHTABLE<int> getconnectedpairs, float xh, float relh, ref int facesplit)
		{
		  // static int timer = NgProfiler::CreateTimer ("AdFront3::GetLocals");
		  // NgProfiler::RegionTimer reg (timer);
        
        
		  if (hashon && faces.Size() < 500)
		  {
			  hashon = 0;
		  }
		  if (hashon && !hashcreated)
		  {
			  hashtable.Create();
			  hashcreated = 1;
		  }
        
		  int i;
		  int j;
		  PointIndex pstind = new PointIndex();
		  Point3d midp = new Point3d();
		  Point3d p0 = new Point3d();
        
		  //  static Array<int, PointIndex::BASE> invpindex;
        
		  Array<MiniElement2d> locfaces2 = new Array<MiniElement2d>(); //all local faces in radius xh
		  Array<int> locfaces3 = new Array<int>(); // all faces in outer radius relh
		  Array<int> findex2 = new Array<int>();
        
		  locfaces2.SetSize(0);
		  locfaces3.SetSize(0);
		  findex2.SetSize(0);
        
		  int cluster = faces.Get(fstind).cluster;
        
		  pstind = faces.Get(fstind).Face().PNum(1);
		  p0 = points[pstind].P();
        
		  locfaces2.Append(faces.Get(fstind).Face());
		  findex2.Append(fstind);
        
        
		  Box3d b1 = new Box3d(p0 - new Vec3d(xh, xh, xh), p0 + new Vec3d(xh, xh, xh));
        
		  if (hashon)
		  {
			  hashtable.GetLocals(locfaces2, findex2, fstind, p0, xh);
		  }
		  else
		  {
			  for (i = 1; i <= faces.Size(); i++)
			  {
			  MiniElement2d face = faces.Get(i).Face();
			  if (faces.Get(i).cluster == cluster && faces.Get(i).Valid() && i != fstind)
			  {
				  Box3d b2 = new Box3d();
				  b2.SetPoint(points[face[0]].P());
				  b2.AddPoint(points[face[1]].P());
				  b2.AddPoint(points[face[2]].P());
        
				  if (b1.Intersect(b2) != 0)
				  {
				  locfaces2.Append(faces.Get(i).Face());
				  findex2.Append(i);
				  }
			  }
			  }
		  }
        
		  //local faces for inner radius:
		  for (i = 1; i <= locfaces2.Size(); i++)
		  {
			  MiniElement2d face = locfaces2.Get(i);
			  Point3d p1 = points[face[0]].P();
			  Point3d p2 = points[face[1]].P();
			  Point3d p3 = points[face[2]].P();
        
			  midp = Center(p1, p2, p3);
        
			  if (Dist2(midp, p0) <= relh * relh || i == 1)
			  {
				  locfaces.Append(locfaces2.Get(i));
			  findex.Append(findex2.Get(i));
			  }
			  else
			  {
			locfaces3.Append(i);
			  }
		  }
        
		  facesplit = locfaces.Size();
        
        
		  //local faces for outer radius:
		  for (i = 1; i <= locfaces3.Size(); i++)
		  {
			  locfaces.Append(locfaces2.Get(locfaces3.Get(i)));
			  findex.Append(findex2.Get(locfaces3.Get(i)));
		  }
        
        
		  invpindex.SetSize(points.Size());
		  for (i = 1; i <= locfaces.Size(); i++)
		  {
			for (j = 1; j <= locfaces.Get(i).GetNP(); j++)
			{
			PointIndex pi = locfaces.Get(i).PNum(j);
			invpindex[pi] = -1;
			}
		  }
        
		  for (i = 1; i <= locfaces.Size(); i++)
		  {
			  for (j = 1; j <= locfaces.Get(i).GetNP(); j++)
			  {
			  PointIndex pi = locfaces.Get(i).PNum(j);
			  if (invpindex[pi] == -1)
			  {
				  pindex.Append(pi);
					  locpoints.Append(points[pi].P());
				  invpindex[pi] = pindex.Size() - 1 + PointIndex.BASE;
			  }
				  // locfaces.Elem(i).PNum(j) = locpoints.Append (points[pi].P());
				  // }
			  // else
				  locfaces.Elem(i).PNum(j) = invpindex[pi];
			  }
		  }
        
        
        
		  if (connectedpairs)
		  {
			  for (i = 1; i <= locpoints.Size(); i++)
			  {
			  int pind = pindex.Get(i);
			  if (pind >= 1 && pind <= connectedpairs.Size())
			  {
				  for (j = 1; j <= connectedpairs.EntrySize(pind); j++)
				  {
				  int oi = connectedpairs.Get(pind, j);
				  int other = invpindex.Get(oi);
				  if (other >= 1 && other <= pindex.Size() && pindex.Get(other) == oi)
				  {
					  // INDEX_2 coned(i, other);
					  // coned.Sort();
					  // (*testout) << "connected: " << locpoints.Get(i) << "-" << locpoints.Get(other) << endl;
					  getconnectedpairs.Set(INDEX_2.Sort(i, other), 1);
				  }
				  }
			  }
			  }
		  }
        
        
		  /*
		    // add isolated points
		  for (i = 1; i <= points.Size(); i++)
		    if (points.Elem(i).Valid() && Dist (points.Elem(i).P(), p0) <= xh)
		      {
			if (!invpindex.Get(i))
			  {
				locpoints.Append (points.Get(i).P());
				pindex.Append (i);
				invpindex.Elem(i) = pindex.Size();
			  }
		      }
		      */
		  return faces.Get(fstind).QualClass();
		}
	}
}