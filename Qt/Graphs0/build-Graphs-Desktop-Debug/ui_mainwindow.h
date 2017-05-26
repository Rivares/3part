/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.2.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QFrame>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTextBrowser>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *action;
    QAction *action_3;
    QWidget *centralWidget;
    QFrame *line;
    QLabel *label;
    QLabel *outputGraph;
    QWidget *layoutWidget;
    QVBoxLayout *verticalLayout_6;
    QVBoxLayout *verticalLayout;
    QVBoxLayout *verticalLayout_5;
    QVBoxLayout *verticalLayout_4;
    QVBoxLayout *verticalLayout_3;
    QHBoxLayout *horizontalLayout_4;
    QHBoxLayout *horizontalLayout;
    QLabel *label_2;
    QLineEdit *inputLeftX;
    QHBoxLayout *horizontalLayout_3;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_3;
    QLineEdit *inputRightX;
    QLabel *label_4;
    QHBoxLayout *horizontalLayout_5;
    QHBoxLayout *horizontalLayout_6;
    QLabel *label_5;
    QLineEdit *inputLeftY;
    QHBoxLayout *horizontalLayout_7;
    QHBoxLayout *horizontalLayout_8;
    QLabel *label_6;
    QLineEdit *inputRightY;
    QLabel *label_7;
    QHBoxLayout *horizontalLayout_9;
    QLabel *label_8;
    QLineEdit *inputAccuracy;
    QFrame *line_3;
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *horizontalLayout_10;
    QRadioButton *Hexch;
    QLabel *label_9;
    QHBoxLayout *horizontalLayout_11;
    QRadioButton *Mexch;
    QLabel *label_10;
    QPushButton *draw;
    QPushButton *clear;
    QPushButton *save;
    QPushButton *exit;
    QTextBrowser *valV1;
    QTextBrowser *valV2;
    QTextBrowser *valV3;
    QFrame *line_4;
    QLabel *label_11;
    QTextBrowser *valF1;
    QTextBrowser *valF2;
    QTextBrowser *valF3;
    QLabel *label_12;
    QLabel *label_13;
    QFrame *line_2;
    QMenuBar *menuBar;
    QMenu *menu;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QStringLiteral("MainWindow"));
        MainWindow->resize(818, 526);
        action = new QAction(MainWindow);
        action->setObjectName(QStringLiteral("action"));
        action_3 = new QAction(MainWindow);
        action_3->setObjectName(QStringLiteral("action_3"));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        line = new QFrame(centralWidget);
        line->setObjectName(QStringLiteral("line"));
        line->setGeometry(QRect(550, 20, 20, 391));
        line->setFrameShape(QFrame::VLine);
        line->setFrameShadow(QFrame::Sunken);
        label = new QLabel(centralWidget);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(640, 0, 81, 21));
        outputGraph = new QLabel(centralWidget);
        outputGraph->setObjectName(QStringLiteral("outputGraph"));
        outputGraph->setGeometry(QRect(20, 30, 531, 371));
        layoutWidget = new QWidget(centralWidget);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(580, 30, 216, 325));
        verticalLayout_6 = new QVBoxLayout(layoutWidget);
        verticalLayout_6->setSpacing(6);
        verticalLayout_6->setContentsMargins(11, 11, 11, 11);
        verticalLayout_6->setObjectName(QStringLiteral("verticalLayout_6"));
        verticalLayout_6->setContentsMargins(0, 0, 0, 0);
        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(6);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        verticalLayout_5 = new QVBoxLayout();
        verticalLayout_5->setSpacing(6);
        verticalLayout_5->setObjectName(QStringLiteral("verticalLayout_5"));
        verticalLayout_4 = new QVBoxLayout();
        verticalLayout_4->setSpacing(6);
        verticalLayout_4->setObjectName(QStringLiteral("verticalLayout_4"));
        verticalLayout_3 = new QVBoxLayout();
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(6);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        label_2 = new QLabel(layoutWidget);
        label_2->setObjectName(QStringLiteral("label_2"));

        horizontalLayout->addWidget(label_2);

        inputLeftX = new QLineEdit(layoutWidget);
        inputLeftX->setObjectName(QStringLiteral("inputLeftX"));
        inputLeftX->setReadOnly(true);

        horizontalLayout->addWidget(inputLeftX);


        horizontalLayout_4->addLayout(horizontalLayout);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        label_3 = new QLabel(layoutWidget);
        label_3->setObjectName(QStringLiteral("label_3"));

        horizontalLayout_2->addWidget(label_3);

        inputRightX = new QLineEdit(layoutWidget);
        inputRightX->setObjectName(QStringLiteral("inputRightX"));

        horizontalLayout_2->addWidget(inputRightX);


        horizontalLayout_3->addLayout(horizontalLayout_2);

        label_4 = new QLabel(layoutWidget);
        label_4->setObjectName(QStringLiteral("label_4"));

        horizontalLayout_3->addWidget(label_4);


        horizontalLayout_4->addLayout(horizontalLayout_3);


        verticalLayout_3->addLayout(horizontalLayout_4);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setSpacing(6);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        label_5 = new QLabel(layoutWidget);
        label_5->setObjectName(QStringLiteral("label_5"));

        horizontalLayout_6->addWidget(label_5);

        inputLeftY = new QLineEdit(layoutWidget);
        inputLeftY->setObjectName(QStringLiteral("inputLeftY"));
        inputLeftY->setEnabled(true);
        inputLeftY->setReadOnly(false);

        horizontalLayout_6->addWidget(inputLeftY);


        horizontalLayout_5->addLayout(horizontalLayout_6);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setSpacing(6);
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setSpacing(6);
        horizontalLayout_8->setObjectName(QStringLiteral("horizontalLayout_8"));
        label_6 = new QLabel(layoutWidget);
        label_6->setObjectName(QStringLiteral("label_6"));

        horizontalLayout_8->addWidget(label_6);

        inputRightY = new QLineEdit(layoutWidget);
        inputRightY->setObjectName(QStringLiteral("inputRightY"));

        horizontalLayout_8->addWidget(inputRightY);


        horizontalLayout_7->addLayout(horizontalLayout_8);

        label_7 = new QLabel(layoutWidget);
        label_7->setObjectName(QStringLiteral("label_7"));

        horizontalLayout_7->addWidget(label_7);


        horizontalLayout_5->addLayout(horizontalLayout_7);


        verticalLayout_3->addLayout(horizontalLayout_5);


        verticalLayout_4->addLayout(verticalLayout_3);

        horizontalLayout_9 = new QHBoxLayout();
        horizontalLayout_9->setSpacing(6);
        horizontalLayout_9->setObjectName(QStringLiteral("horizontalLayout_9"));
        label_8 = new QLabel(layoutWidget);
        label_8->setObjectName(QStringLiteral("label_8"));

        horizontalLayout_9->addWidget(label_8);

        inputAccuracy = new QLineEdit(layoutWidget);
        inputAccuracy->setObjectName(QStringLiteral("inputAccuracy"));

        horizontalLayout_9->addWidget(inputAccuracy);


        verticalLayout_4->addLayout(horizontalLayout_9);


        verticalLayout_5->addLayout(verticalLayout_4);

        line_3 = new QFrame(layoutWidget);
        line_3->setObjectName(QStringLiteral("line_3"));
        line_3->setFrameShape(QFrame::HLine);
        line_3->setFrameShadow(QFrame::Sunken);

        verticalLayout_5->addWidget(line_3);

        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        horizontalLayout_10 = new QHBoxLayout();
        horizontalLayout_10->setSpacing(6);
        horizontalLayout_10->setObjectName(QStringLiteral("horizontalLayout_10"));
        Hexch = new QRadioButton(layoutWidget);
        Hexch->setObjectName(QStringLiteral("Hexch"));
        Hexch->setMinimumSize(QSize(55, 20));
        Hexch->setMaximumSize(QSize(55, 16777215));
        QFont font;
        font.setItalic(false);
        font.setUnderline(false);
        font.setKerning(true);
        Hexch->setFont(font);
        Hexch->setChecked(true);
        Hexch->setAutoRepeat(false);

        horizontalLayout_10->addWidget(Hexch);

        label_9 = new QLabel(layoutWidget);
        label_9->setObjectName(QStringLiteral("label_9"));

        horizontalLayout_10->addWidget(label_9);


        verticalLayout_2->addLayout(horizontalLayout_10);

        horizontalLayout_11 = new QHBoxLayout();
        horizontalLayout_11->setSpacing(6);
        horizontalLayout_11->setObjectName(QStringLiteral("horizontalLayout_11"));
        Mexch = new QRadioButton(layoutWidget);
        Mexch->setObjectName(QStringLiteral("Mexch"));
        Mexch->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_11->addWidget(Mexch);

        label_10 = new QLabel(layoutWidget);
        label_10->setObjectName(QStringLiteral("label_10"));

        horizontalLayout_11->addWidget(label_10);


        verticalLayout_2->addLayout(horizontalLayout_11);


        verticalLayout_5->addLayout(verticalLayout_2);


        verticalLayout->addLayout(verticalLayout_5);

        draw = new QPushButton(layoutWidget);
        draw->setObjectName(QStringLiteral("draw"));

        verticalLayout->addWidget(draw);

        clear = new QPushButton(layoutWidget);
        clear->setObjectName(QStringLiteral("clear"));

        verticalLayout->addWidget(clear);

        save = new QPushButton(layoutWidget);
        save->setObjectName(QStringLiteral("save"));

        verticalLayout->addWidget(save);

        exit = new QPushButton(layoutWidget);
        exit->setObjectName(QStringLiteral("exit"));

        verticalLayout->addWidget(exit);


        verticalLayout_6->addLayout(verticalLayout);

        valV1 = new QTextBrowser(centralWidget);
        valV1->setObjectName(QStringLiteral("valV1"));
        valV1->setGeometry(QRect(220, 420, 71, 21));
        valV2 = new QTextBrowser(centralWidget);
        valV2->setObjectName(QStringLiteral("valV2"));
        valV2->setGeometry(QRect(330, 420, 71, 21));
        valV3 = new QTextBrowser(centralWidget);
        valV3->setObjectName(QStringLiteral("valV3"));
        valV3->setGeometry(QRect(450, 420, 61, 21));
        line_4 = new QFrame(centralWidget);
        line_4->setObjectName(QStringLiteral("line_4"));
        line_4->setGeometry(QRect(10, 400, 551, 20));
        line_4->setFrameShape(QFrame::HLine);
        line_4->setFrameShadow(QFrame::Sunken);
        label_11 = new QLabel(centralWidget);
        label_11->setObjectName(QStringLiteral("label_11"));
        label_11->setGeometry(QRect(40, 440, 91, 20));
        valF1 = new QTextBrowser(centralWidget);
        valF1->setObjectName(QStringLiteral("valF1"));
        valF1->setGeometry(QRect(220, 450, 71, 21));
        valF2 = new QTextBrowser(centralWidget);
        valF2->setObjectName(QStringLiteral("valF2"));
        valF2->setGeometry(QRect(330, 450, 71, 21));
        valF3 = new QTextBrowser(centralWidget);
        valF3->setObjectName(QStringLiteral("valF3"));
        valF3->setGeometry(QRect(450, 450, 61, 21));
        label_12 = new QLabel(centralWidget);
        label_12->setObjectName(QStringLiteral("label_12"));
        label_12->setGeometry(QRect(170, 420, 31, 20));
        label_12->setStyleSheet(QStringLiteral(""));
        label_13 = new QLabel(centralWidget);
        label_13->setObjectName(QStringLiteral("label_13"));
        label_13->setGeometry(QRect(150, 450, 61, 20));
        line_2 = new QFrame(centralWidget);
        line_2->setObjectName(QStringLiteral("line_2"));
        line_2->setGeometry(QRect(560, 20, 212, 3));
        line_2->setFrameShape(QFrame::HLine);
        line_2->setFrameShadow(QFrame::Sunken);
        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 818, 27));
        menu = new QMenu(menuBar);
        menu->setObjectName(QStringLiteral("menu"));
        MainWindow->setMenuBar(menuBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        MainWindow->setStatusBar(statusBar);

        menuBar->addAction(menu->menuAction());
        menu->addAction(action);
        menu->addSeparator();
        menu->addAction(action_3);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", 0));
#ifndef QT_NO_TOOLTIP
        MainWindow->setToolTip(QApplication::translate("MainWindow", "<html><head/><body><p>Simulation Program</p></body></html>", 0));
#endif // QT_NO_TOOLTIP
        action->setText(QApplication::translate("MainWindow", "\320\241\320\276\321\205\321\200\320\260\320\275\320\270\321\202\321\214 \320\272\320\260\321\200\321\202\320\270\320\275\320\272\321\203", 0));
        action_3->setText(QApplication::translate("MainWindow", "\320\222\321\213\321\205\320\276\320\264", 0));
        label->setText(QApplication::translate("MainWindow", "\320\230\320\275\321\202\320\265\321\200\320\262\320\260\320\273\321\213:", 0));
        outputGraph->setText(QApplication::translate("MainWindow", "TextLabel", 0));
        label_2->setText(QApplication::translate("MainWindow", "t = [", 0));
        inputLeftX->setText(QApplication::translate("MainWindow", "0", 0));
        label_3->setText(QApplication::translate("MainWindow", ",", 0));
        inputRightX->setText(QApplication::translate("MainWindow", "1000", 0));
        label_4->setText(QApplication::translate("MainWindow", "]", 0));
        label_5->setText(QApplication::translate("MainWindow", "T 'C = [", 0));
        inputLeftY->setText(QApplication::translate("MainWindow", "100", 0));
        label_6->setText(QApplication::translate("MainWindow", ",", 0));
        inputRightY->setText(QApplication::translate("MainWindow", "180", 0));
        label_7->setText(QApplication::translate("MainWindow", "]", 0));
        label_8->setText(QApplication::translate("MainWindow", "\320\242\320\276\321\207\320\275\320\276\321\201\321\202\321\214", 0));
        inputAccuracy->setText(QApplication::translate("MainWindow", "700", 0));
        Hexch->setText(QString());
        label_9->setText(QApplication::translate("MainWindow", "<html><head/><body><p>\320\242\320\265\320\277\320\273\320\276\320\276\320\261\320\274\320\265\320\275\320\275\320\260\321\217 \320\274\320\276\320\264\320\265\320\273\321\214</p></body></html>", 0));
        Mexch->setText(QString());
        label_10->setText(QApplication::translate("MainWindow", "<html><head/><body><p>\320\234\320\260\321\201\321\201\320\276\320\276\320\261\320\274\320\265\320\275\320\275\320\260\321\217 \320\274\320\276\320\264\320\265\320\273\321\214</p></body></html>", 0));
        draw->setText(QApplication::translate("MainWindow", "\320\237\320\276\321\201\321\202\321\200\320\276\320\270\321\202\321\214", 0));
        clear->setText(QApplication::translate("MainWindow", "\320\236\321\202\321\207\320\270\321\201\321\202\320\270\321\202\321\214", 0));
        save->setText(QApplication::translate("MainWindow", "\320\241\320\276\321\205\321\200\320\260\320\275\320\270\321\202\321\214", 0));
        exit->setText(QApplication::translate("MainWindow", "\320\222\321\213\320\271\321\202\320\270", 0));
        label_11->setText(QApplication::translate("MainWindow", "\320\243\321\201\321\202\320\260\320\275\320\276\320\262.\320\267\320\275\320\260\321\207.", 0));
        label_12->setText(QApplication::translate("MainWindow", "<html><head/><body><p><span style=\" color:#299d32;\">\320\237\320\260\321\200:</span></p></body></html>", 0));
        label_13->setText(QApplication::translate("MainWindow", "<html><head/><body><p><span style=\" color:#2677cd;\">\320\226\320\270\320\264\320\272\320\276\321\201\321\202\321\214:</span></p></body></html>", 0));
        menu->setTitle(QApplication::translate("MainWindow", "\320\244\320\260\320\271\320\273", 0));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
