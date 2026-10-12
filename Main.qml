import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.15
import factorizer 1.0

ApplicationWindow {
    id: window
    width: 800
    height: 650
    visible: true
    title: "Параллельный факторизатор целых чисел"

    // Экземпляр бэкенда из C++
    FactorizerBackend {
        id: backend
    }

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 20
        spacing: 14

        RowLayout {
            Layout.fillWidth: true
            spacing: 15

            Label {
                text: "Разрядность (бит):"
                font.bold: true
                Layout.alignment: Qt.AlignVCenter
            }

            ComboBox {
                id: bitCombo
                model: [128, 256, 512, 1024, 2048, 4096]
                currentIndex: 0
                enabled: !backend.isRunning
            }

            Button {
                text: "Случайное число"
                enabled: !backend.isRunning
                onClicked: {
                    inputField.text = backend.generateRandomNumber(bitCombo.currentValue)
                }
            }
        }

        Label {
            text: "Число для факторизации:"
            font.bold: true
        }

        ScrollView {
            Layout.fillWidth: true
            Layout.preferredHeight: 100

            TextArea {
                id: inputField
                placeholderText: "Введите число..."
                wrapMode: TextEdit.WrapAnywhere
                font.family: "Monospace"
                text: "172694645192493905907229620525860152663"
            }
        }

        RowLayout {
            Layout.fillWidth: true
            spacing: 12

            Button {
                text: "Факторизовать"
                highlighted: true
                enabled: !backend.isRunning && inputField.text.trim().length > 0
                onClicked: {
                    backend.startFactorization(inputField.text.trim(), bitCombo.currentValue)
                }
            }

            Button {
                text: "Остановить"
                enabled: backend.isRunning
                onClicked: backend.stop()
            }

            BusyIndicator {
                running: backend.isRunning
                Layout.preferredHeight: 34
                Layout.preferredWidth: 34
            }

            Label {
                text: backend.statusText
                font.italic: true
                Layout.alignment: Qt.AlignVCenter
            }
        }

        Label {
            text: "Результат:"
            font.bold: true
        }

        ScrollView {
            Layout.fillWidth: true
            Layout.fillHeight: true

            TextArea {
                id: resultArea
                text: backend.resultText
                readOnly: true
                wrapMode: TextEdit.WrapAnywhere
                font.family: "Monospace"
            }
        }
    }
}