import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.15
import factorizer 1.0

ApplicationWindow {
    id: window

    // Переключатель темы оформления (по умолчанию тёмная)
    // Свойство напрямую берется из бэкенда:
    property bool darkMode: backend.darkMode
    // Динамическая палитра цветов
    readonly property color colBg: darkMode ? "#18191c" : "#f4f6f9"
    readonly property color colCardBg: darkMode ? "#24262b" : "#ffffff"
    readonly property color colBorder: darkMode ? "#383c45" : "#cfd8dc"
    readonly property color colBorderFocus: darkMode ? "#4da3ff" : "#1976d2"
    readonly property color colTextPrimary: darkMode ? "#f0f2f5" : "#212529"
    readonly property color colTextSecondary: darkMode ? "#9ba1ad" : "#546e7a"
    readonly property color colBadgeBg: darkMode ? "#1f334d" : "#e3f2fd"
    readonly property color colBadgeBorder: darkMode ? "#2b5687" : "#90caf9"
    readonly property color colBadgeText: darkMode ? "#70b5ff" : "#1565c0"
    readonly property int inputBitLength: backend.getBitLength(sanitizeNumber(inputField.text))

    function sanitizeNumber(str) {
        return str.replace(/\s+/g, "").trim();
    }

    function formatDigitsByGroups(str, groupSize) {
        var clean = sanitizeNumber(str);
        if (clean.length === 0)
            return "";

        var size = groupSize || 10;
        var formatted = "";
        for (var i = 0; i < clean.length; i += size) {
            if (i > 0)
                formatted += " ";

            formatted += clean.substring(i, Math.min(i + size, clean.length));
        }
        return formatted;
    }

    // Именованная функция автоподбора разрядности
    function updateBitSelection() {
        if (!autoBitCheckBox.checked || inputBitLength <= 0 || backend.isRunning)
            return ;

        var bitModels = bitCombo.model;
        for (var i = 0; i < bitModels.length; ++i) {
            if (inputBitLength <= bitModels[i]) {
                if (bitCombo.currentIndex !== i)
                    bitCombo.currentIndex = i;

                return ;
            }
        }
        bitCombo.currentIndex = bitModels.length - 1;
    }

    width: 520
    height: 840
    minimumWidth: 420
    minimumHeight: 680
    visible: true
    title: "Параллельный факторизатор целых чисел"
    color: colBg
    onInputBitLengthChanged: {
        updateBitSelection();
    }

    FactorizerBackend {
        id: backend

        onIsRunningChanged: {
            if (backend.isRunning) {
                elapsedTimer.startTime = Date.now();
                elapsedTimer.elapsedSeconds = 0;
                elapsedTimer.running = true;
            } else {
                elapsedTimer.running = false;
            }
        }
    }

    Timer {
        id: elapsedTimer

        property var startTime: 0
        property double elapsedSeconds: 0

        function formatElapsed() {
            var totalSec = Math.floor(elapsedSeconds);
            var mins = Math.floor(totalSec / 60);
            var secs = totalSec % 60;
            var tenths = Math.floor((elapsedSeconds - totalSec) * 10);
            var minsStr = (mins < 10 ? "0" : "") + mins;
            var secsStr = (secs < 10 ? "0" : "") + secs;
            return minsStr + ":" + secsStr + "." + tenths;
        }

        interval: 100
        repeat: true
        running: false
        onTriggered: {
            elapsedSeconds = (Date.now() - startTime) / 1000;
        }
    }

    // Вспомогательный элемент для надёжного копирования любого текста в буфер обмена
    TextEdit {
        id: clipboardHelper

        function copyText(val) {
            text = val;
            selectAll();
            copy();
        }

        visible: false
    }

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 16
        spacing: 12

        // --- 1. Панель параметров темы ---
        RowLayout {
            Layout.fillWidth: true
            spacing: 8

            Label {
                text: "Параметры факторизатора"
                font.bold: true
                font.pixelSize: 15
                color: colTextPrimary
                Layout.alignment: Qt.AlignVCenter
            }

            Item {
                Layout.fillWidth: true
            }

            Switch {
                id: themeSwitch

                text: backend.darkMode ? "🌙 Тёмная" : "☀️ Светлая"
                checked: backend.darkMode
                onToggled: backend.darkMode = checked

                contentItem: Text {
                    text: themeSwitch.text
                    font: themeSwitch.font
                    color: colTextPrimary
                    verticalAlignment: Text.AlignVCenter
                    leftPadding: themeSwitch.indicator.width + themeSwitch.spacing
                }

            }

        }

        // --- 2. Разрядность (без кнопки случайного числа) ---
        RowLayout {
            Layout.fillWidth: true
            spacing: 8

            Label {
                text: "Разрядность:"
                color: colTextPrimary
                font.pixelSize: 13
                Layout.alignment: Qt.AlignVCenter
            }

            ComboBox {
                id: bitCombo

                model: [128, 256, 512, 1024, 2048, 4096]
                currentIndex: 0
                implicitWidth: 90
                enabled: !backend.isRunning && !autoBitCheckBox.checked
            }

            CheckBox {
                id: autoBitCheckBox

                text: "Авто"
                checked: true
                enabled: !backend.isRunning
                Layout.alignment: Qt.AlignVCenter
                onCheckedChanged: {
                    if (checked)
                        updateBitSelection();

                }

                contentItem: Text {
                    text: autoBitCheckBox.text
                    font: autoBitCheckBox.font
                    color: autoBitCheckBox.enabled ? colTextPrimary : colTextSecondary
                    verticalAlignment: Text.AlignVCenter
                    leftPadding: autoBitCheckBox.indicator.width + autoBitCheckBox.spacing
                }

            }

            Item {
                Layout.fillWidth: true
            }

        }

        // --- 3. Заголовок ввода и кнопки быстрого действия ---
        RowLayout {
            Layout.fillWidth: true
            spacing: 8

            Label {
                text: "Входное число:"
                font.bold: true
                color: colTextPrimary
            }

            Label {
                text: {
                    if (window.inputBitLength > 0) {
                        var digits = sanitizeNumber(inputField.text).replace(/^0+/, "").length;
                        return "(" + window.inputBitLength + " бит, " + digits + " цифр)";
                    }
                    return "(0 бит)";
                }
                font.pixelSize: 12
                font.bold: true
                color: window.inputBitLength > bitCombo.currentValue ? "#e53935" : (window.darkMode ? "#66bb6a" : "#2e7d32")
                Layout.alignment: Qt.AlignVCenter
            }

            Item {
                Layout.fillWidth: true
            }

            Button {
                text: "📋 Вставить"
                enabled: !backend.isRunning
                implicitHeight: 28
                onClicked: {
                    inputField.selectAll();
                    inputField.paste();
                    if (groupDigitsSwitch.checked)
                        inputField.text = formatDigitsByGroups(inputField.text, 10);

                }
            }

            Button {
                text: "✕ Стереть"
                enabled: !backend.isRunning && inputField.text.length > 0
                implicitHeight: 28
                onClicked: inputField.clear()
            }

        }

        // --- 4. Контрастное поле ввода ---
        ScrollView {
            Layout.fillWidth: true
            Layout.preferredHeight: 120

            TextArea {
                id: inputField

                placeholderText: "Введите десятичное целое число..."
                wrapMode: TextEdit.WrapAnywhere
                font.family: "Monospace"
                font.pixelSize: 13
                selectByMouse: true
                color: colTextPrimary
                leftPadding: 10
                rightPadding: 10
                topPadding: 8
                bottomPadding: 8
                text: "172694645192493905907229620525860152663"

                background: Rectangle {
                    color: colCardBg
                    radius: 6
                    border.width: inputField.activeFocus ? 2 : 1
                    border.color: inputField.activeFocus ? colBorderFocus : colBorder
                }

            }

        }

        // Переключатель группировки цифр по 10
        RowLayout {
            Layout.fillWidth: true

            Switch {
                id: groupDigitsSwitch

                text: "Группировать цифры по 10"
                checked: false
                enabled: !backend.isRunning
                onToggled: {
                    var raw = sanitizeNumber(inputField.text);
                    if (checked)
                        inputField.text = formatDigitsByGroups(raw, 10);
                    else
                        inputField.text = raw;
                }

                contentItem: Text {
                    text: groupDigitsSwitch.text
                    font: groupDigitsSwitch.font
                    color: groupDigitsSwitch.enabled ? colTextPrimary : colTextSecondary
                    verticalAlignment: Text.AlignVCenter
                    leftPadding: groupDigitsSwitch.indicator.width + groupDigitsSwitch.spacing
                }

            }

            Item {
                Layout.fillWidth: true
            }

        }

        // --- 5. Панель действий: Случайное число, Запуск/Остановка, Таймер ---
        RowLayout {
            Layout.fillWidth: true
            spacing: 8

            // Кнопка генерации случайного числа
            Button {
                text: "🎲 Случайное число"
                implicitHeight: 36
                enabled: !backend.isRunning
                onClicked: {
                    var rnd = backend.generateRandomNumber(bitCombo.currentValue);
                    if (groupDigitsSwitch.checked)
                        inputField.text = formatDigitsByGroups(rnd, 10);
                    else
                        inputField.text = rnd;
                }
            }

            // Динамическая кнопка запуска/остановки
            Button {
                id: actionButton

                Layout.fillWidth: true
                implicitHeight: 36
                text: backend.isRunning ? "⏹ Остановить" : "▶ Факторизовать"
                highlighted: !backend.isRunning
                enabled: backend.isRunning || (sanitizeNumber(inputField.text).length > 0 && window.inputBitLength <= 4096)
                onClicked: {
                    if (backend.isRunning)
                        backend.stop();
                    else
                        backend.startFactorization(sanitizeNumber(inputField.text), bitCombo.currentValue);
                }

                background: Rectangle {
                    radius: 6
                    color: {
                        if (!actionButton.enabled)
                            return window.darkMode ? "#2b2d33" : "#e0e0e0";

                        if (backend.isRunning)
                            return actionButton.down ? "#b71c1c" : (window.darkMode ? "#c62828" : "#d32f2f");

                        return actionButton.down ? "#1565c0" : (window.darkMode ? "#1976d2" : "#1976d2");
                    }
                }

                contentItem: Text {
                    text: actionButton.text
                    font.bold: true
                    font.pixelSize: 13
                    color: actionButton.enabled ? "#ffffff" : (window.darkMode ? "#616161" : "#9e9e9e")
                    horizontalAlignment: Text.AlignHCenter
                    verticalAlignment: Text.AlignVCenter
                }

            }

            BusyIndicator {
                running: backend.isRunning
                visible: backend.isRunning
                Layout.preferredHeight: 28
                Layout.preferredWidth: 28
            }

            Rectangle {
                visible: elapsedTimer.running || elapsedTimer.elapsedSeconds > 0
                implicitWidth: timerLabel.implicitWidth + 14
                implicitHeight: 28
                radius: 4
                color: colBadgeBg
                border.color: colBadgeBorder
                Layout.alignment: Qt.AlignVCenter

                Label {
                    id: timerLabel

                    anchors.centerIn: parent
                    text: (backend.isRunning ? "⏳ " : "⏱ ") + elapsedTimer.formatElapsed() + " с"
                    font.family: "Monospace"
                    font.bold: true
                    font.pixelSize: 11
                    color: colBadgeText
                }

            }

        }

        // Строка статуса текущего шага
        RowLayout {
            Layout.fillWidth: true

            Label {
                text: backend.statusText.length > 0 ? backend.statusText : "Готов к работе"
                font.italic: true
                color: colTextSecondary
                font.pixelSize: 12
                Layout.fillWidth: true
                elide: Text.ElideRight
                verticalAlignment: Text.AlignVCenter
            }

        }

        // --- 6. Заголовок результата с кнопкой Копировать ---
        RowLayout {
            Layout.fillWidth: true
            spacing: 8

            Label {
                text: "Результат факторизации:"
                font.bold: true
                color: colTextPrimary
            }

            Item {
                Layout.fillWidth: true
            }

            Button {
                id: copyBtn

                text: "📋 Копировать"
                enabled: backend.resultText.length > 0
                implicitHeight: 28
                onClicked: {
                    clipboardHelper.copyText(backend.resultText);
                    copyAnim.restart();
                }

                SequentialAnimation {
                    id: copyAnim

                    PropertyAction {
                        target: copyBtn
                        property: "text"
                        value: "✓ Скопировано!"
                    }

                    PauseAnimation {
                        duration: 1500
                    }

                    PropertyAction {
                        target: copyBtn
                        property: "text"
                        value: "📋 Копировать"
                    }

                }

            }

        }

        // --- 7. Контрастное поле вывода результата ---
        ScrollView {
            Layout.fillWidth: true
            Layout.fillHeight: true

            TextArea {
                id: resultArea

                text: backend.resultText
                readOnly: true
                wrapMode: TextEdit.WrapAnywhere
                font.family: "Monospace"
                font.pixelSize: 12
                color: colTextPrimary
                selectByMouse: true
                leftPadding: 10
                rightPadding: 10
                topPadding: 8
                bottomPadding: 8

                background: Rectangle {
                    color: colCardBg
                    radius: 6
                    border.width: 1
                    border.color: colBorder
                }

            }

        }

    }

}
