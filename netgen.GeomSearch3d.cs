namespace netgen
{

	public class GeomSearch3d
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void AddElem(MiniElement2d elem, int elemnum)
		  {
			Point3d minp = new Point3d();
			Point3d maxp = new Point3d();
			ElemMaxExt(minp, maxp, elem);
			int sx = (int)(minp.X() - minext.X()) / elemsize.X() + 1.0;
			int ex = (int)(maxp.X() - minext.X()) / elemsize.X() + 1.0;
			int sy = (int)(minp.Y() - minext.Y()) / elemsize.Y() + 1.0;
			int ey = (int)(maxp.Y() - minext.Y()) / elemsize.Y() + 1.0;
			int sz = (int)(minp.Z() - minext.Z()) / elemsize.Z() + 1.0;
			int ez = (int)(maxp.Z() - minext.Z()) / elemsize.Z() + 1.0;
        
			for (int ix = sx; ix <= ex; ix++)
			{
			  for (int iy = sy; iy <= ey; iy++)
			  {
				for (int iz = sz; iz <= ez; iz++)
				{
					int ind = ix + (iy - 1) * size.i1 + (iz - 1) * size.i2 * size.i1;
					if (ind < 1 || ind > size.i1 * size.i2 * size.i3)
					{
						cerr << "Illegal hash-position";
						cerr << "Position: " << ix << "," << iy << "," << iz << "\n";
					throw new Exception("Illegal position in Geomsearch");
					}
					hashtable.Elem(ind).Append(elemnum);
				}
			  }
			}
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GetLocals(Array<MiniElement2d> locfaces, Array<int> findex, int fstind, Point3d p0, double xh)
		  {
			hashcount++;
        
			Point3d minp = new Point3d();
			Point3d maxp = new Point3d();
			Point3d midp = new Point3d();
        
			minp.CopyFrom(p0 - new Vec3d(xh, xh, xh)); //lay cube over sphere
			maxp.CopyFrom(p0 + new Vec3d(xh, xh, xh));
        
			MaxCoords(minext,minp); //cube may not be out of hash-region
			MinCoords(maxextreal,maxp);
        
        
			int cluster = faces.Get(fstind).Cluster();
        
			int sx = (int)(minp.X() - minext.X()) / elemsize.X() + 1.0;
			int ex = (int)(maxp.X() - minext.X()) / elemsize.X() + 1.0;
			int sy = (int)(minp.Y() - minext.Y()) / elemsize.Y() + 1.0;
			int ey = (int)(maxp.Y() - minext.Y()) / elemsize.Y() + 1.0;
			int sz = (int)(minp.Z() - minext.Z()) / elemsize.Z() + 1.0;
			int ez = (int)(maxp.Z() - minext.Z()) / elemsize.Z() + 1.0;
			int ix;
			int iy;
			int iz;
			int i;
			int k;
        
			int cnt1 = 0; // test, how efficient hashtable is
			int cnt2 = 0;
			int cnt3 = 0;
        
			for (ix = sx; ix <= ex; ix++)
			{
			for (iy = sy; iy <= ey; iy++)
			{
				for (iz = sz; iz <= ez; iz++)
				{
				int ind = ix + (iy - 1) * size.i1 + (iz - 1) * size.i2 * size.i1;
        
				//go through all elements in one hash area
				Array<int> area = *hashtable.Elem(ind);
				for (k = 1; k <= area.Size(); k++)
				{
					cnt2++;
					i = area.Get(k);
					if (faces.Get(i).Cluster() == cluster && faces.Get(i).Valid() && faces.Get(i).HashValue() != hashcount && i != fstind)
					{
					cnt1++;
					MiniElement2d face = faces.Get(i).Face();
        
					Point3d p1 = (*points)[face.PNum(1)].P();
					Point3d p2 = (*points)[face.PNum(2)].P();
					Point3d p3 = (*points)[face.PNum(3)].P();
        
					midp = Center(p1, p2, p3);
        
					// if (Dist2 (midp, p0) <= xh*xh)  
								if ((Dist2(p1, p0) <= xh * xh) || (Dist2(p2, p0) <= xh * xh) || (Dist2(p3, p0) <= xh * xh) || (Dist2(midp, p0) <= xh * xh)) // by Jochen Wild
								{
						cnt3++;
						locfaces.Append(faces.Get(i).Face());
						findex.Append(i);
						faces.Elem(i).SetHashValue(hashcount);
								}
					}
				}
				}
			}
			}
			/*
			  if (faces->Size() != 0 && hashcount % 200 == 0)
			  {
			  (*mycout) << "n.o.f= " << faces->Size();
			  (*mycout) << ", n.o.lf= " << locfaces.Size();
			  (*mycout) << ", hashf= " << (double)cnt2/(double)faces->Size();
			  (*mycout) << " (" << (double)cnt1/(double)faces->Size();
			  (*mycout) << ", " << (double)cnt3/(double)faces->Size() << ")" << endl;
			  }
			*/
        
		  }
	}
}