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
    readonly property color colWarning: darkMode ? "#ffb74d" : "#b26a00"
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

        // Ближайшая разрядность, в которую помещается число (64..512)
        bitSpin.value = Math.min(Math.max(inputBitLength, bitSpin.from), bitSpin.to);
    }

    // 700: при 520 ряды кнопок и переключателей (включая выбор метода)
    // не помещались и обрезались справа.
    width: 700
    height: 840
    minimumWidth: 560
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

            SpinBox {
                id: bitSpin

                // Произвольная разрядность 64..512. Ограничения по разрядности
                // у методов сняты: при явном выборе SIQS/GNFS расчёт запустится
                // и на 512 битах (пользователь сам решает, что готов ждать).
                from: 64
                to: 512
                stepSize: 1
                value: 128
                editable: true
                implicitWidth: 115
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

            // Генерация сложного числа p×q (semiprime) заданной разрядности
            Button {
                text: "🔐 Полупростое число: p×q"
                implicitHeight: 28
                enabled: !backend.isRunning
                Layout.alignment: Qt.AlignVCenter
                onClicked: backend.generateSemiprime(bitSpin.value)
            }

        }

        // Получение сгенерированного сложного числа из бэкенда
        Connections {
            target: backend

            function onSemiprimeReady(number) {
                if (groupDigitsSwitch.checked)
                    inputField.text = formatDigitsByGroups(number, 10);
                else
                    inputField.text = number;
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
                color: window.inputBitLength > bitSpin.value ? "#e53935" : (window.darkMode ? "#66bb6a" : "#2e7d32")
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

            // Переключатель пропуска алгоритма ECM
            Switch {
                id: skipEcmSwitch

                text: "Пропускать ECM"
                checked: backend.skipEcm
                enabled: !backend.isRunning
                onToggled: backend.skipEcm = checked

                contentItem: Text {
                    text: skipEcmSwitch.text
                    font: skipEcmSwitch.font
                    color: skipEcmSwitch.enabled ? colTextPrimary : colTextSecondary
                    verticalAlignment: Text.AlignVCenter
                    leftPadding: skipEcmSwitch.indicator.width + skipEcmSwitch.spacing
                }

                ToolTip.visible: hovered
                ToolTip.text: "Полностью отключает лестницу ECM (для полупростых p×q:\n"
                              + "сразу работают решёта SIQS/GNFS). Методы Полларда — работают.\n\n"
                              + "ECM применяется только в режиме «Авто» к числам от 136 бит:\n"
                              + "меньшие числа мгновенно раскладывает быстрый SIQS.\n"
                              + "При явном выборе метода (SIQS/GNFS) ECM не запускается и так."
            }

            // Явный выбор метода решета: авто / SIQS / GNFS
            Label {
                text: "Метод:"
                color: colTextPrimary
                font.pixelSize: 13
                Layout.alignment: Qt.AlignVCenter
            }

            ComboBox {
                id: methodCombo

                model: ["Авто", "SIQS", "GNFS"]
                currentIndex: backend.method
                implicitWidth: 118
                implicitHeight: 28
                font.pixelSize: 12
                enabled: !backend.isRunning
                Layout.alignment: Qt.AlignVCenter
                onActivated: backend.method = currentIndex

                contentItem: Text {
                    leftPadding: 8
                    rightPadding: methodCombo.indicator.width + methodCombo.spacing
                    text: methodCombo.displayText
                    font: methodCombo.font
                    color: methodCombo.enabled ? colTextPrimary : colTextSecondary
                    verticalAlignment: Text.AlignVCenter
                    elide: Text.ElideRight
                }

                background: Rectangle {
                    radius: 6
                    color: window.darkMode ? "#2b2d33" : "#ffffff"
                    border.width: 1
                    border.color: methodCombo.enabled
                                  ? (methodCombo.activeFocus ? colBorderFocus : colBorder)
                                  : colBorder
                }

                ToolTip.visible: hovered
                ToolTip.text: "Метод факторизации:\n"
                              + "• Авто — SIQS до 320 бит, GNFS в полосе 128-150 бит;\n"
                              + "• SIQS — всегда квадратичное решето, без ограничения разрядности;\n"
                              + "• GNFS — всегда решето числового поля, без ограничения разрядности.\n\n"
                              + "При явном выборе работает только выбранный метод (плюс пробное\n"
                              + "деление, Поллард и проверка простоты): лестница ECM не запускается,\n"
                              + "фолбэка на другой метод нет.\n\n"
                              + "Запустить 512-битное число можно любым методом, но SIQS выше ~350 бит\n"
                              + "практически не набирает соотношений (программа сразу выведет оценку\n"
                              + "сложности и предупреждение), а GNFS этой реализации рассчитан\n"
                              + "примерно до 150 бит."
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
                    var rnd = backend.generateRandomNumber(bitSpin.value);
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
                enabled: backend.isRunning || (sanitizeNumber(inputField.text).length > 0 && window.inputBitLength <= bitSpin.to)
                onClicked: {
                    if (backend.isRunning)
                        backend.stop();
                    else
                        backend.startFactorization(sanitizeNumber(inputField.text), bitSpin.value);
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

        // Постоянное предупреждение о неприменимости метода: строка статуса
        // обновляется прогрессом, а эта — нет (см. warningText в бэкенде).
        Label {
            Layout.fillWidth: true
            visible: backend.warningText.length > 0
            text: backend.warningText
            color: colWarning
            font.pixelSize: 12
            wrapMode: Text.WordWrap
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
