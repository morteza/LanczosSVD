#include <QtTest/QtTest>

class TestEVD: public QObject
{
    Q_OBJECT
private slots:
    void doTest();
};

void TestEVD::doTest()
{
    //QCOMPARE(QString("TEST"), QString("TEST"));
}

//QTEST_MAIN(TestEVD)
#include "TestEVD.moc"
