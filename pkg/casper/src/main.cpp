// qasper2c.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "casper.h"

void example(Casper* c)
{	
	//Exon* e1 = new Exon(1, 100, 400);
	//Exon* e2 = new Exon(2, 1000, 1100);
	//Exon* e3 = new Exon(3, 2000, 2500);
	Exon* e1 = new Exon(1, 300);
	Exon* e2 = new Exon(2, 100);
	Exon* e3 = new Exon(3, 500);
	c->addExon(e1);
	c->addExon(e2);
	c->addExon(e3);

	Variant* v1 = new Variant(1, 3);
	v1->addExon(e1);
	v1->addExon(e2);
	v1->addExon(e3);
	Variant* v2 = new Variant(2, 2);
	v2->addExon(e1);
	v2->addExon(e3);
	Variant* v3 = new Variant(3, 2);
	v3->addExon(e1);
	v3->addExon(e2);

Fragment* f0 = new Fragment(1, 1, 4069);
f0->left[0] = 3;
f0->right[0] = 3;
c->addFragment(f0);
Fragment* f1 = new Fragment(1, 1, 292);
f1->left[0] = 1;
f1->right[0] = 3;
c->addFragment(f1);
Fragment* f2 = new Fragment(2, 1, 654);
f2->left[0] = 2;
f2->left[1] = 3;
f2->right[0] = 3;
c->addFragment(f2);
Fragment* f3 = new Fragment(1, 1, 1807);
f3->left[0] = 1;
f3->right[0] = 1;
c->addFragment(f3);
Fragment* f4 = new Fragment(1, 2, 974);
f4->left[0] = 1;
f4->right[0] = 1;
f4->right[1] = 2;
c->addFragment(f4);
Fragment* f5 = new Fragment(2, 2, 399);
f5->left[0] = 1;
f5->left[1] = 2;
f5->right[0] = 2;
f5->right[1] = 3;
c->addFragment(f5);
Fragment* f6 = new Fragment(2, 1, 246);
f6->left[0] = 1;
f6->left[1] = 2;
f6->right[0] = 3;
c->addFragment(f6);
Fragment* f7 = new Fragment(1, 1, 227);
f7->left[0] = 2;
f7->right[0] = 3;
c->addFragment(f7);
Fragment* f8 = new Fragment(2, 1, 366);
f8->left[0] = 1;
f8->left[1] = 3;
f8->right[0] = 3;
c->addFragment(f8);
Fragment* f9 = new Fragment(1, 2, 225);
f9->left[0] = 1;
f9->right[0] = 2;
f9->right[1] = 3;
c->addFragment(f9);
Fragment* f10 = new Fragment(1, 2, 396);
f10->left[0] = 1;
f10->right[0] = 1;
f10->right[1] = 3;
c->addFragment(f10);
Fragment* f11 = new Fragment(1, 1, 339);
f11->left[0] = 1;
f11->right[0] = 2;
c->addFragment(f11);
Fragment* f12 = new Fragment(2, 1, 3);
f12->left[0] = 1;
f12->left[1] = 2;
f12->right[0] = 2;
c->addFragment(f12);
Fragment* f13 = new Fragment(2, 2, 1);
f13->left[0] = 2;
f13->left[1] = 3;
f13->right[0] = 2;
f13->right[1] = 3;
c->addFragment(f13);
Fragment* f14 = new Fragment(2, 2, 1);
f14->left[0] = 1;
f14->left[1] = 2;
f14->right[0] = 1;
f14->right[1] = 2;
c->addFragment(f14);
Fragment* f15 = new Fragment(1, 2, 1);
f15->left[0] = 2;
f15->right[0] = 2;
f15->right[1] = 3;
c->addFragment(f15);


	c->addVariant(v1);
	c->addVariant(v2);
	c->addVariant(v3);
}
void txample(Casper* c)
{
	FILE *fs;
	char line[1024];

	fs = fopen("exons.dat", "rt");
	while (fgets(line, 1024, fs) != NULL)
	{
		int eid;
		int length;
		sscanf(line, "%d\t%d", &eid, &length);
		Exon* e = new Exon(eid, length);
		c->addExon(e);
	}
	fclose(fs);

	fs = fopen("fragments.dat", "rt");
	while (fgets(line, 1024, fs) != NULL)
	{
		int leftc;
		int rightc;
		int count;
		sscanf(line, "%d\t%d\t%d", &leftc, &rightc, &count);
		Fragment* f = new Fragment(leftc, rightc, count);
		for (int l = 0; l < leftc; l++)
		{
			int eid;
			fgets(line, 1024, fs);
			sscanf(line, "%d", &eid);
			f->left[l] = eid;
		}
		for (int r = 0; r < rightc; r++)
		{
			int eid;
			fgets(line, 1024, fs);
			sscanf(line, "%d", &eid);
			f->right[r] = eid;
		}
		c->addFragment(f);
	}
	fclose(fs);

	fs = fopen("variants.dat", "rt");
	while (fgets(line, 1024, fs) != NULL)
	{
		int id;
		int exonc;
		sscanf(line, "%d\t%d", &id, &exonc);
		Variant* v = new Variant(id, exonc);
		for (int l = 0; l < exonc; l++)
		{
			int eid;
			fgets(line, 1024, fs);
			sscanf(line, "%d", &eid);
			v->addExon(c->getExon(eid));
		}
		c->addVariant(v);
	}
	fclose(fs);
}

DiscreteDF* readlenDist()
{
	DiscreteDF* ld = new DiscreteDF(300);

	ld->Set(101, 0.000117373600169339);
	ld->Set(102, 8.42682257626024e-05);
	ld->Set(103, 0.000101322509547891);
	ld->Set(104, 6.82171351411543e-05);
	ld->Set(105, 7.12267146326758e-05);
	ld->Set(106, 9.1290577909486e-05);
	ld->Set(107, 8.22618394349214e-05);
	ld->Set(108, 6.21979761581113e-05);
	ld->Set(109, 7.82490667795593e-05);
	ld->Set(110, 5.01596581920252e-05);
	ld->Set(111, 7.02235214688353e-05);
	ld->Set(112, 7.12267146326758e-05);
	ld->Set(113, 7.92522599433998e-05);
	ld->Set(114, 6.01915898304303e-05);
	ld->Set(115, 5.41724308473872e-05);
	ld->Set(116, 7.72458736157188e-05);
	ld->Set(117, 4.61468855366632e-05);
	ld->Set(118, 5.81852035027493e-05);
	ld->Set(119, 6.92203283049948e-05);
	ld->Set(120, 3.91245333897797e-05);
	ld->Set(121, 4.71500787005037e-05);
	ld->Set(122, 6.11947829942708e-05);
	ld->Set(123, 6.52075556496328e-05);
	ld->Set(124, 3.61149538982582e-05);
	ld->Set(125, 5.41724308473872e-05);
	ld->Set(126, 4.41404992089822e-05);
	ld->Set(127, 5.21660445197062e-05);
	ld->Set(128, 4.71500787005037e-05);
	ld->Set(129, 3.61149538982582e-05);
	ld->Set(130, 4.51436923728227e-05);
	ld->Set(131, 4.51436923728227e-05);
	ld->Set(132, 6.82171351411543e-05);
	ld->Set(133, 6.01915898304303e-05);
	ld->Set(134, 5.81852035027493e-05);
	ld->Set(135, 6.32011693219518e-05);
	ld->Set(136, 6.82171351411543e-05);
	ld->Set(137, 8.92841915818049e-05);
	ld->Set(138, 9.32969642371669e-05);
	ld->Set(139, 7.02235214688353e-05);
	ld->Set(140, 8.12586462710809e-05);
	ld->Set(141, 9.02873847456454e-05);
	ld->Set(142, 0.000101322509547891);
	ld->Set(143, 0.000107341668530934);
	ld->Set(144, 0.000141450236101511);
	ld->Set(145, 0.000146466201920714);
	ld->Set(146, 0.000150478974576076);
	ld->Set(147, 0.000157501326722959);
	ld->Set(148, 0.000181577962655131);
	ld->Set(149, 0.000225718461864113);
	ld->Set(150, 0.000230734427683316);
	ld->Set(151, 0.000268855767909255);
	ld->Set(152, 0.000299954755988311);
	ld->Set(153, 0.000276881313219979);
	ld->Set(154, 0.00035813995949106);
	ld->Set(155, 0.000511628513558657);
	ld->Set(156, 0.000494574229773369);
	ld->Set(157, 0.000568810523897566);
	ld->Set(158, 0.000710260759999077);
	ld->Set(159, 0.000690196896722267);
	ld->Set(160, 0.000789513019942477);
	ld->Set(161, 0.00081459284903849);
	ld->Set(162, 0.000920931324405583);
	ld->Set(163, 0.00111354441186296);
	ld->Set(164, 0.00136133312333156);
	ld->Set(165, 0.00153990150649517);
	ld->Set(166, 0.00176662316152313);
	ld->Set(167, 0.00221505050575983);
	ld->Set(168, 0.00248290308050525);
	ld->Set(169, 0.00288819311869681);
	ld->Set(170, 0.00349211540332880);
	ld->Set(171, 0.0039676289629892);
	ld->Set(172, 0.0045795767929319);
	ld->Set(173, 0.00534300679061453);
	ld->Set(174, 0.00628500517146076);
	ld->Set(175, 0.00693507434162941);
	ld->Set(176, 0.00813589655874649);
	ld->Set(177, 0.00905181191733287);
	ld->Set(178, 0.00988145266382897);
	ld->Set(179, 0.0109960002688558);
	ld->Set(180, 0.0127014286473846);
	ld->Set(181, 0.0131458432189660);
	ld->Set(182, 0.0145723838979472);
	ld->Set(183, 0.0157812316603750);
	ld->Set(184, 0.0162868410149506);
	ld->Set(185, 0.0174846536525761);
	ld->Set(186, 0.0179280650309937);
	ld->Set(187, 0.0186483577226311);
	ld->Set(188, 0.019328522687715);
	ld->Set(189, 0.0197980170883924);
	ld->Set(190, 0.0200628600836462);
	ld->Set(191, 0.0207480410145493);
	ld->Set(192, 0.0219368249137003);
	ld->Set(193, 0.0206386929596907);
	ld->Set(194, 0.0220000260830223);
	ld->Set(195, 0.0215385572276556);
	ld->Set(196, 0.0209657339311027);
	ld->Set(197, 0.0222648690782762);
	ld->Set(198, 0.0222849329415530);
	ld->Set(199, 0.0212386024716673);
	ld->Set(200, 0.0207380090829109);
	ld->Set(201, 0.0207701112641538);
	ld->Set(202, 0.0202705210685612);
	ld->Set(203, 0.0200999782307083);
	ld->Set(204, 0.0198822853141550);
	ld->Set(205, 0.0187847919929134);
	ld->Set(206, 0.0191950979969242);
	ld->Set(207, 0.0183012528879423);
	ld->Set(208, 0.0177595285794684);
	ld->Set(209, 0.0178728894069824);
	ld->Set(210, 0.0170693316827462);
	ld->Set(211, 0.016321952775685);
	ld->Set(212, 0.0172749862813335);
	ld->Set(213, 0.0162748026969845);
	ld->Set(214, 0.0155204014377764);
	ld->Set(215, 0.0152806382716186);
	ld->Set(216, 0.0143988314806028);
	ld->Set(217, 0.0137076313907167);
	ld->Set(218, 0.0136975994590782);
	ld->Set(219, 0.0129833259264238);
	ld->Set(220, 0.0127816841004919);
	ld->Set(221, 0.0122560108826394);
	ld->Set(222, 0.0117514047212277);
	ld->Set(223, 0.0111274185733189);
	ld->Set(224, 0.0109859683372174);
	ld->Set(225, 0.0101844169993088);
	ld->Set(226, 0.00966977890625862);
	ld->Set(227, 0.00936179860495959);
	ld->Set(228, 0.00904077679253062);
	ld->Set(229, 0.00821615201185373);
	ld->Set(230, 0.00767743728287138);
	ld->Set(231, 0.00711263953162918);
	ld->Set(232, 0.00683274863891768);
	ld->Set(233, 0.00612549745841012);
	ld->Set(234, 0.00567205414835421);
	ld->Set(235, 0.00520155655451301);
	ld->Set(236, 0.00476015156242319);
	ld->Set(237, 0.00408700894948621);
	ld->Set(238, 0.00364961673005175);
	ld->Set(239, 0.00320720854479809);
	ld->Set(240, 0.00283803346050479);
	ld->Set(241, 0.0024247178770025);
	ld->Set(242, 0.00219498664248302);
	ld->Set(243, 0.00199234162338724);
	ld->Set(244, 0.00156297494926351);
	ld->Set(245, 0.00135431077118468);
	ld->Set(246, 0.00120583818293629);
	ld->Set(247, 0.000951027119320798);
	ld->Set(248, 0.000885819563671165);
	ld->Set(249, 0.000712267146326758);
	ld->Set(250, 0.000584861614519014);
	ld->Set(251, 0.000486548684462645);
	ld->Set(252, 0.000415321969829969);
	ld->Set(253, 0.00039726449288084);
	ld->Set(254, 0.000286913244858384);
	ld->Set(255, 0.000308983494462875);
	ld->Set(256, 0.000249795097796286);
	ld->Set(257, 0.000171546031016726);
	ld->Set(258, 0.000150478974576076);
	ld->Set(259, 0.000168536451525205);
	ld->Set(260, 0.000140447042937671);
	ld->Set(261, 0.000111354441186296);
	ld->Set(262, 0.000121386372824701);
	ld->Set(263, 0.000117373600169339);
	ld->Set(264, 8.52714189264429e-05);
	ld->Set(265, 8.42682257626024e-05);
	ld->Set(266, 0.000116370407005499);
	ld->Set(267, 9.7309736892529e-05);
	ld->Set(268, 6.01915898304303e-05);
	ld->Set(269, 7.42362941241973e-05);
	ld->Set(270, 7.92522599433998e-05);
	ld->Set(271, 8.02554531072404e-05);
	ld->Set(272, 8.62746120902834e-05);
	ld->Set(273, 6.92203283049948e-05);
	ld->Set(274, 7.02235214688353e-05);
	ld->Set(275, 6.11947829942708e-05);
	ld->Set(276, 6.72139419773138e-05);
	ld->Set(277, 6.21979761581113e-05);
	ld->Set(278, 6.01915898304303e-05);
	ld->Set(279, 5.71820103389087e-05);
	ld->Set(280, 5.41724308473872e-05);
	ld->Set(281, 5.01596581920252e-05);
	ld->Set(282, 5.61788171750682e-05);
	ld->Set(283, 4.41404992089822e-05);
	ld->Set(284, 2.80894085875341e-05);
	ld->Set(285, 5.01596581920252e-05);
	ld->Set(286, 5.11628513558657e-05);
	ld->Set(287, 4.91564650281847e-05);
	ld->Set(288, 4.41404992089822e-05);
	ld->Set(289, 3.91245333897797e-05);
	ld->Set(290, 3.31053744067366e-05);
	ld->Set(291, 3.61149538982582e-05);
	ld->Set(292, 2.80894085875341e-05);
	ld->Set(293, 3.71181470620987e-05);
	ld->Set(294, 4.31373060451417e-05);
	ld->Set(295, 2.50798290960126e-05);
	ld->Set(296, 2.90926017513746e-05);
	ld->Set(297, 4.31373060451417e-05);
	ld->Set(298, 2.80894085875341e-05);
	ld->Set(299, 3.71181470620987e-05);
	return ld;
}

double cumu_fragsta(double x)
{
	return x;
}
int main()
{
	Casper* c = new Casper();
	c->setDistributions(readlenDist(), cumu_fragsta);
	txample(c);

	double* pi = c->calculateEM(c->randomPi());
	for (int i = 0; i < c->varcount(); i++)
	{
		printf("%f\n", pi[i]);
	}
	getchar();

	return 0;
}
