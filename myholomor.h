#pragma once
#include <igl/cotmatrix.h>
#include "mos.h"
#include <Eigen\QR>
//#include <Eigen\SparseQR>

//#include "shortest.h"
class Oneform
{
public:
	Oneform(PMyMesh &mesh1, Mos &mos1) :mesh(mesh1), mos(mos1)
	{
		printf_s("create Oneform \r\n");
	}
	~Oneform()
	{}
	void init_g2()
	{
		g_2 = 2 + mesh.n_edges() - mesh.n_vertices() - mesh.n_faces();
		printf_s("g_2= %d\r\n", g_2);
	}
	void init_w()
	{
		for (PMyMesh::FaceIter fter = mesh.faces_begin(); fter != mesh.faces_end(); fter++)
		{
			mesh.data(*fter).w.resize(3, g_2);
			mesh.data(*fter).w_star.resize(3, g_2);
			mesh.data(*fter).w.setZero();
			mesh.data(*fter).w_star.setZero();
		}
	}
	float weight(EdgeHandle e)
	{
		HalfedgeHandle he = mesh.s_halfedge_handle(e, 0);
		he = mesh.next_halfedge_handle(he);
		double angle1 = mesh.calc_sector_angle(he);
		he = mesh.s_halfedge_handle(e, 1);
		he = mesh.next_halfedge_handle(he);
		double angle2 = mesh.calc_sector_angle(he);
		double temp = 1 / 2.0*(cos(angle1) / sin(angle1) + cos(angle2) / sin(angle2));
		//temp = angle1;
		mesh.data(e).weight = temp;
		return temp;
	}
	void init_weight()
	{
		igl::cotmatrix(V, F, L);
		HalfedgeHandle he;
		for (PMyMesh::EdgeIter eter = mesh.edges_begin(); eter != mesh.edges_end(); eter++)
		{
			he = mesh.s_halfedge_handle(*eter, 0);
			mesh.data(*eter).weight = L.coeff(mesh.to_vertex_handle(he).idx(), mesh.from_vertex_handle(he).idx());
			printf_s("%f,  %f\r\n", mesh.data(*eter).weight, weight(*eter));
			if (L.coeff(mesh.to_vertex_handle(he).idx(), mesh.from_vertex_handle(he).idx()) != L.coeff(mesh.from_vertex_handle(he).idx(), mesh.to_vertex_handle(he).idx()))
			{
				printf_s("错误\r\n");
			}
		}
	}
	void one_omiga(VertexHandle v, int i)
	{
		g_2 = 2 + mesh.n_edges() - mesh.n_vertices() - mesh.n_faces();
		HalfedgeHandle he = mesh.halfedge_handle(v), he1;
		EdgeHandle e = mesh.edge_handle(he);
		VertexHandle  v1 = v;
		//PMyMesh::Point p;
		while (!mesh.data(e).sharp_)
		{
			he = mesh.ccw_rotated_halfedge_handle(he);
			e = mesh.edge_handle(he);
		}
		if (!mesh.data(e).sharp_)
		{
			printf_s("错误\r\n");
			return;
		}
		if (mesh.from_vertex_handle(he) != v)
		{
			printf_s("错误\r\n");
			return;
		}
		int times = 0;
		do {
			times++;
			printf_s("times: %d\r\n", times);

			he = mesh.ccw_rotated_halfedge_handle(he);
			//openmesh中关于判断是否he为null的等价判断
			while (!mesh.data(mesh.edge_handle(he)).sharp_&&mesh.is_valid_handle(mesh.face_handle(he))) {

				mesh.data(he).omiga[i] = -1;
				mesh.data(mesh.opposite_halfedge_handle(he)).omiga[i] = 1;
				mesh.data(mesh.edge_handle(he)).sharp1_ = true;
				he = mesh.ccw_rotated_halfedge_handle(he);
				//e = mesh.edge_handle(he);
			}
			if (!mesh.data(mesh.edge_handle(he)).sharp_)
			{
				printf_s("错误\r\n");
				return;
			}
			he = mesh.opposite_halfedge_handle(he);
		} while (mesh.from_vertex_handle(he) != v1);
	}
	void set_omiga()
	{
		for (int i = 0; i < g_2; i++)
		{
			mos.init_sharp();
			mos.findsaddle(mos.vh[i]);
			one_omiga(mos.vh[i], i);
		}
		/*mos.init_A();
		mos.set_f();
		mos.setkey();*/

	}
	void harmonic_oneform(int ii)
	{

		Eigen::SparseMatrix<double> A(mesh.n_vertices(), mesh.n_vertices());
		A.setZero();
		Eigen::VectorXd bx(mesh.n_vertices()), x(mesh.n_vertices());
		bx.setZero();
		//by.setZero();
		HalfedgeHandle he; int i = 0, j = 0;
		for (PMyMesh::EdgeIter eter = mesh.edges_begin(); eter != mesh.edges_end(); eter++)
		{
			he = mesh.s_halfedge_handle(*eter, 0);
			i = mesh.to_vertex_handle(he).idx();
			j = mesh.from_vertex_handle(he).idx();
			A.coeffRef(i, j) = mesh.data(*eter).weight;
			A.coeffRef(j, i) = mesh.data(*eter).weight;
			A.coeffRef(i, i) -= mesh.data(*eter).weight;
			A.coeffRef(j, j) -= mesh.data(*eter).weight;
			bx.coeffRef(i) -= mesh.data(*eter).weight*mesh.data(mesh.opposite_halfedge_handle(he)).omiga[ii];
			bx.coeffRef(j) -= mesh.data(*eter).weight*mesh.data(he).omiga[ii];
			//printf_s("%f  \r\n", mesh.data(*eter).weight*mesh.data(mesh.opposite_halfedge_handle(he)).omiga[ii]);

			if (mesh.data(mesh.opposite_halfedge_handle(he)).omiga[ii] != -mesh.data(he).omiga[ii])
			{
				//printf_s("cuowu\r\n");
				return;
			}
			//printf_s("A: %f, %f,b: %f, %f, %f, %f\r\n", A.coeff(i, i), A.coeff(j, j), bx.coeff(i), bx.coeff(j), mesh.data(*eter).weight*mesh.data(mesh.opposite_halfedge_handle(he)).omiga[ii], mesh.data(*eter).weight*mesh.data(he).omiga[ii]);
		}
		for (int Ti = 0; Ti < bx.rows(); Ti++)
		{
			printf_s("bx: %e\r\n", bx.coeff(Ti));
		}
		A = A * 1;
		bx = bx * 1;
		for (int Ti = 0; Ti < bx.rows(); Ti++)
		{
			//printf_s("new bx: %e\r\n", bx.coeff(Ti));
		}

		A.makeCompressed();
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>ldlt;

		printf_s("compressed make\r\n");
		ldlt.compute(A);
		x = ldlt.solve(bx);
		Eigen::VectorXd temp = A * x;
		for (int Ti = 0; Ti < temp.rows(); Ti++)
		{
			//printf_s("x: %f, temp: %e,bx: %e\r\n", x.coeff(Ti), temp.coeff(Ti), bx.coeff(Ti));
		}
		//printf_s("%d  %d\r\n",temp.rows(),mesh.n_vertices());
		for (PMyMesh::HalfedgeIter hter = mesh.halfedges_begin(); hter != mesh.halfedges_end(); hter++)
		{
			i = mesh.to_vertex_handle(*hter).idx();
			//printf_s("%d\r\n",i);
			j = mesh.from_vertex_handle(*hter).idx();
			mesh.data(*hter).omiga[ii] += (x.coeff(i) - x.coeff(j));
			//printf_s("x %f\r\n",x.coeff(i));
		}
	}
	void set_harmonic_oneform()
	{
		for (int i = 0; i < g_2; i++)
		{
			harmonic_oneform(i);
		}
	}
	void derivate_basis(FaceHandle f, int ii)
	{
		HalfedgeHandle he[3];
		he[0] = mesh.halfedge_handle(f);
		he[1] = mesh.next_halfedge_handle(he[0]);
		he[2] = mesh.prev_halfedge_handle(he[0]);
		PMyMesh::Point p[3];
		for (int i = 0; i < 3; i++)
		{
			p[i] = mesh.point(mesh.to_vertex_handle(he[i]));
		}
		Eigen::Matrix3f m3;
		Eigen::Vector3f n;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++) { m3.coeffRef(j, i) = p[i][j]; }
		}
		Vec3f nn = mesh.calc_face_normal(f);
		n.coeffRef(0) = nn[0];
		n.coeffRef(1) = nn[1];
		n.coeffRef(2) = nn[2];
		//n = ((m3.col(1) - m3.col(0)).cross(m3.col(2)-m3.col(1)));
		//n.normalize();
		//printf_s("n:%f,%f,%f , normal: %f,%f,%f\r\n",n(0),n(1),n(2), (mesh.calc_face_normal(f))[0], (mesh.calc_face_normal(f))[1], (mesh.calc_face_normal(f))[2]);

		Eigen::Vector3f temp;
		for (int i = 0; i < 3; i++)
		{
			temp.coeffRef(i) = mesh.data(he[(i + 2) % 3]).omiga[ii];
		}
		//printf_s("temp: %f,n.norm: %f, n*v: %f \r\n",temp(0)+temp(1)+temp(2),n.norm(),(n.transpose()*(m3.col(2) - m3.col(1)))(0));
		temp = m3 * temp;
		mesh.data(f).w.col(ii) = temp.cross(n);
		mesh.data(f).w.col(ii) /= (2.0* mesh.calc_sector_area(he[0]));
		printf_s("omiga %f, test: %f\r\n", mesh.data(he[2]).omiga[ii], ((m3.col(2) - m3.col(1)).transpose()*mesh.data(f).w.col(ii)).coeff(0));
		//printf_s("test: %f\r\n", (n.transpose()*mesh.data(f).w.col(ii))(0));
		//return -1 * temp.cross(n);
	}
	void one_cotangentbasis(int ii)
	{
		for (PMyMesh::FaceIter fter = mesh.faces_begin(); fter != mesh.faces_end(); fter++)
		{
			derivate_basis(*fter, ii);
		}
	}
	void set_cotangentbasis()
	{
		for (int i = 0; i < g_2; i++)
		{
			one_cotangentbasis(i);
		}
	}
	float inner(Eigen::MatrixXf &w1, Eigen::MatrixXf &w2)
	{
		float sum = 0; int id = 0;
		printf_s("w1,%d,%d,w2: %d,%d\r\n", w1.cols(), w1.rows(), w2.cols(), w2.rows());
		Eigen::Vector3f n, temp1, temp2; Vec3f nn; float area;
		for (PMyMesh::FaceIter fter = mesh.faces_begin(); fter != mesh.faces_end(); fter++)
		{
			
			id = (*fter).idx();
			area = mesh.calc_sector_area(mesh.halfedge_handle(*fter));
			if (area < 0)
			{
				printf_s("cuowu\r\n");
				return 0;
			}
			temp1 = (w1.row(id)).transpose();
			temp2 = (w2.row(id)).transpose();
			sum += (temp1.transpose()*temp2).coeff(0)/area;
			//sum += ((temp1.transpose()*temp2).coeff(0)+(temp1.transpose()*mesh.data(*fter).w.col(2)).coeff(0)+ (temp2.transpose()*mesh.data(*fter).w.col(2)).coeff(0)+ mesh.data(*fter).w.col(2).norm());
			//printf_s("test: %f\r\n", (n.transpose()*(temp2)).coeff(0));
			//sum+=(((w1.row(id).cross(w.row(id)))*temp)*mesh.calc_sector_area(mesh.halfedge_handle(*fter))).coeff(0);
			//sum+=(((w1.row(id)).cross(w2.row(id)))*(temp*mesh.calc_sector_area(mesh.halfedge_handle(*fter)))).coeff(0);
		}
		return sum;
	}
	void one_cotangentstar(int ii)
	{
		Eigen::Vector3f n, temp; Vec3f nn;
		for (PMyMesh::FaceIter fter = mesh.faces_begin(); fter != mesh.faces_end(); fter++)
		{
			nn = mesh.calc_face_normal(*fter);
			n.coeffRef(0) = nn[0];
			n.coeffRef(1) = nn[1];
			n.coeffRef(2) = nn[2];
			temp = mesh.data(*fter).w.col(ii);
			mesh.data(*fter).w_star.col(ii) = n.cross(temp);
			printf_s("w_star norm: %f, w norm:: %f\r\n", mesh.data(*fter).w_star.col(ii).norm(), mesh.data(*fter).w.col(ii).norm());
		}
	}
	void set_cotangentstar()
	{
		for (int i = 0; i < g_2; i++)
		{
			one_cotangentstar(i);
		}
	}
	void one_tromigastar(int ii)
	{
		HalfedgeHandle he; VertexHandle v1, v2; Eigen::RowVector3f temp;
		FaceHandle f; float sum = 0; //PMyMesh::Point p;
		for (PMyMesh::EdgeIter eter = mesh.edges_begin(); eter != mesh.edges_end(); eter++)
		{
			sum = 0;
			he = mesh.s_halfedge_handle(*eter, 0);
			v1 = mesh.to_vertex_handle(he);
			v2 = mesh.from_vertex_handle(he);
			temp.coeffRef(0) = mesh.point(v1)[0] - mesh.point(v2)[0];
			temp.coeffRef(1) = mesh.point(v1)[1] - mesh.point(v2)[1];
			temp.coeffRef(2) = mesh.point(v1)[2] - mesh.point(v2)[2];
			f = mesh.face_handle(he);
			sum += (temp*mesh.data(f).w_star.col(ii)).coeff(0);
			f = mesh.face_handle(mesh.opposite_halfedge_handle(he));
			sum += (temp*mesh.data(f).w_star.col(ii)).coeff(0);
			mesh.data(he).omiga_star[ii] = 0.5*sum;
			mesh.data(mesh.opposite_halfedge_handle(he)).omiga_star[ii] = -0.5*sum;
		}
	}
	void one_omigastar(int ii)
	{
		Eigen::VectorXf  b(g_2), x(g_2);
		b.setZero();
		Eigen::MatrixXf w2(mesh.n_faces(), 3), w3(mesh.n_faces(), 3);
		std::vector<Eigen::MatrixXf> w1;
		w2.setZero();
		for (int i = 0; i<g_2; i++)
		{
			Eigen::MatrixXf tempm;
			tempm.resize(mesh.n_faces(), 3);
			for (PMyMesh::FaceIter fter = mesh.faces_begin(); fter != mesh.faces_end(); fter++)
			{
				tempm.row((*fter).idx()) = (mesh.data(*fter).w.col(i)).transpose();
			}
			w1.push_back(tempm);
		}
		for (PMyMesh::FaceIter fter = mesh.faces_begin(); fter != mesh.faces_end(); fter++)
		{
			w2.row((*fter).idx()) = (mesh.data(*fter).w_star.col(ii)).transpose();
		}
		/*for (int i=0;i<w2.rows();i++)
		{
		printf_s("test:%f,test:%f\r\n",(w2.row(i)).norm(),((w1[ii].row(i))).norm());
		}*/
		/*for (int i=0;i<mesh.n_faces();i++)
		{
		printf_s("w[2]:%f,%f,%f,w2:%f,%f,%f\r\n",w1[0].coeff(i,0), w1[0].coeff(i, 1), w1[0].coeff(i, 2),w2.coeff(i,0), w2.coeff(i, 1), w2.coeff(i, 2));
		}*/
		for (int i = 0; i<g_2; i++)
		{
			b.coeffRef(i) = inner(w1[i], w2);
			printf_s("b:%f\r\n", b(i));
		}
		x = ww.inverse()*b;
		for (int i = 0; i<g_2; i++)
		{

			printf("x:%f\r\n",x(i));
		}
		/*w3.setZero();
		for (int i = 0; i<g_2; i++)
		{

			w3 += w1[i] * x(i);
		}
		printf_s("w3:%f\r\n", inner(w2, w3));*/
		//x = ww.lu().solve(b);
		float sum = 0;
		
		for (PMyMesh::EdgeIter eter = mesh.edges_begin(); eter != mesh.edges_end(); eter++)
		{
			sum = 0;
			for (int i = 0; i<g_2; i++)
			{
				sum += mesh.data(mesh.s_halfedge_handle(*eter,0)).omiga[i]*x(i);

			}
			mesh.data(mesh.s_halfedge_handle(*eter, 0)).omiga_star[ii] = sum;
			mesh.data(mesh.s_halfedge_handle(*eter, 1)).omiga_star[ii] = -sum;
		}
	}
	void set_omigastar()
	{
		for (int i = 0; i<g_2; i++)
		{
			one_omigastar(i);
		}

	}

	void init_ww()
	{  
		ww.resize(g_2, g_2);
		ww.setZero();
		
		std::vector<Eigen::MatrixXf> w1;
		for (int i = 0; i<g_2; i++)
		{
			Eigen::MatrixXf tempm;
			tempm.resize(mesh.n_faces(), 3);
			for (PMyMesh::FaceIter fter = mesh.faces_begin(); fter != mesh.faces_end(); fter++)
			{
				tempm.row((*fter).idx()) = (mesh.data(*fter).w.col(i)).transpose();
			}
			w1.push_back(tempm);
		}
		for (int i = 0; i<g_2; i++)
		{
			for (int j = 0; j<g_2; j++)
			{
				
				if (j>=i)
				{
					ww.coeffRef(i, j) = inner(w1[i], w1[j]);
					ww.coeffRef(j, i) = ww.coeff(i, j);
				}
			}
		}
		for (int i = 0; i<g_2; i++)
		{
			for (int j = 0; j<g_2; j++)
			{
				printf_s("ww:%f\r\n", ww.coeff(i, j));
			}
		}
		
	}


	void set_uv(Eigen::MatrixXd &V_uv)
	{
		V_uv.resize(mesh.n_vertices(), 2);
		V_uv.setZero();
		for (PMyMesh::VertexIter vter = mesh.vertices_begin(); vter != mesh.vertices_end(); vter++)
		{
			//V_uv.coeffRef((*vter).idx(), 0) = fabs(sin(mesh.data(vter).v)) / 3.0;
		V_uv.coeffRef((*vter).idx(), 0) = mesh.data(vter).u * 10;
			//V_uv.coeffRef((*vter).idx(), 1) = mesh.data(vter).v;
			//V_uv.coeffRef((*vter).idx(), 1) = fabs(sin(mesh.data(vter).u))/3.0;
			V_uv.coeffRef((*vter).idx(), 1) = mesh.data(vter).v * 10;
			printf_s("u: %f\r\n", mesh.data(vter).v);
		}

	}

	Eigen::MatrixXd V;
	Eigen::SparseMatrix<float> L;
	Eigen::MatrixXi F;
	int g_2;
protected:
	PMyMesh & mesh;

	Mos &mos;
	Eigen::MatrixXf ww;

};
